"""Quadratic-tower field code with eight split charts."""
from functools import lru_cache
from math import comb

M=500
HALF=250
MASK=(1<<M)-1
SUBMASK=(1<<HALF)-1
MOD_LOW=1|(1<<32)|(1<<40)|(1<<59)
MODULUS=(1<<HALF)|MOD_LOW
XI=1<<191
THETA=1<<HALF
LAMBDA=2
THIRD=MASK//3
RADIUS=103
BYTES=64
CHARTS=tuple((outer,j,k) for outer in (0,1) for j in (0,1) for k in (0,1))
FALLBACK=1|(2<<HALF)
ARRAYS=(
    "Y0","Y1","Y2","Y3","KEY0","KEY1","V0","V1","V2","V3","BEST0","BEST1",
    "X2","C1","C2","D0","M0","M1","M2","M3","M4","RATIO","INVCOEF",
    "BASEPIVOT","DELTA","F1","F2","F3","LIN","QUAD",
)
BASE={name:48+8*i for i,name in enumerate(ARRAYS)}
NEXT=48+8*len(ARRAYS)
LAST=NEXT+RADIUS+1
INEXT=LAST+RADIUS+1
ILAST=INEXT+RADIUS+1
MULT=ILAST+RADIUS+1
REDSHIFT1=MULT+16*4
REDSHIFT4=REDSHIFT1+2
XI_MAP=REDSHIFT4+16
MAP_BASE=XI_MAP+32*256*4
MAP_TAGS=(1,"lambda","theta","theta2")

def remainder(a, b=MODULUS):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def carryless(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def base_product(a, b):
    return remainder(carryless(a, b))


def base_power(a, e):
    out = 1
    for bit in bin(e)[2:]:
        out = base_product(out, out)
        if bit == "1":
            out = base_product(out, a)
    return out


def parts(a):
    return a & SUBMASK, a >> HALF


def join(a, b):
    return a | (b << HALF)


def product(a, b):
    a0, a1 = parts(a)
    b0, b1 = parts(b)
    z0 = base_product(a0, b0)
    z1 = base_product(a1, b1)
    z2 = base_product(a0 ^ a1, b0 ^ b1)
    return join(z0 ^ base_product(XI, z1), z2 ^ z0)


def square(a):
    a0, a1 = parts(a)
    z0, z1 = base_product(a0, a0), base_product(a1, a1)
    return join(z0 ^ base_product(XI, z1), z1)


def good(a):
    a0, a1 = parts(a)
    return a0 != 0 and a1 != 0 and a0 != a1


def packed_words(a):
    for part in parts(a):
        for j in range(4):
            yield (part >> (64*j)) & ((1 << 64)-1)



def cube(a):
    return product(a, square(a))


def frobenius(a, s):
    for _ in range(s):
        a = square(a)
    return a


def power(a, e):
    result = 1
    for bit in bin(e)[2:]:
        result = square(result)
        if bit == "1":
            result = product(result, a)
    return result


def codeword(a,b):
    return a,b,product(a,b),cube(a)^product(LAMBDA,cube(b))


def encode(message):
    return sum(x<<(M*i) for i,x in enumerate(codeword(message&MASK,message>>M)))


@lru_cache(None)
def constant(tag):
    return {"lambda":LAMBDA,"theta":THETA,"theta2":square(THETA)}[tag]


def map_columns(tag):
    for half in (0,1):
        for i in range(256):
            if i>=HALF:
                yield 0
                continue
            x=1<<(half*HALF+i)
            yield square(x) if tag==1 else product(constant(tag),x)


def map_words(tag):
    cols = iter(map_columns(tag))
    for _ in range(BYTES):
        vals = [0]
        for _ in range(8):
            v = next(cols)
            vals.extend(x ^ v for x in vals[:])
        for value in vals:
            yield from packed_words(value)


def xi_entry(j, value):
    a = base_product(XI, remainder(value << (8*j)))
    return tuple((a >> (64*k)) & ((1 << 64)-1) for k in range(4))


def chart_row(outer, known_other, known_product):
    lo, hi = parts(outer)
    matrix = ((lo, base_product(XI, hi)), (hi, lo ^ hi))
    return matrix[known_product][known_other], matrix[known_product][1-known_other]


def split_recover(outer, other_index, product_index, p, q):
    t, v = chart_row(outer, other_index, product_index)
    s = base_power(v, SUBMASK-1)
    z = base_product(s, q ^ base_product(t, p))
    other = join(p,z) if other_index == 0 else join(z,p)
    return other, product(outer, other)


def directions(outer, other_index, product_index):
    t, v = chart_row(outer, other_index, product_index)
    s = base_power(v, SUBMASK-1)
    coefficient = base_product(s,t)
    for lane in (0, 1):
        for i in range(HALF):
            e = 1 << i
            p = e if lane == 0 else 0
            z = base_product(coefficient if lane == 0 else s, e)
            other = join(p,z) if other_index == 0 else join(z,p)
            yield other, product(outer,other)


def derivative_columns(alpha, e):
    lin = product(alpha, square(e))
    quad = product(alpha, e)
    for half in (0, 1):
        l, q = (lin, quad) if half == 0 else (
            product(THETA, lin), product(square(THETA), quad))
        for _ in range(HALF):
            yield l ^ q
            l = product(l, 2)
            q = product(q, 4)
        yield from (0 for _ in range(256-HALF))


def derivative_entry(alpha,e,j,value):
    cols = tuple(derivative_columns(alpha,e))
    out = product(alpha,cube(e)) if j == 0 else 0
    for bit in range(8):
        if (value >> bit) & 1:
            out ^= cols[8*j+bit]
    return out


def reduction_table(s):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << s))


DIRECTIONS=MAP_BASE+len(MAP_TAGS)*BYTES*256*8
SAVED=DIRECTIONS+M*16
AFFINE=SAVED+(RADIUS+1)*24
AFFINE_STRIDE=BYTES*256*8
MEMORY_END=AFFINE+M*AFFINE_STRIDE


def W_bound():
    nodes=sum(comb(2*M,i) for i in range(RADIUS+1))
    outer=sum(comb(M,i) for i in range(RADIUS+1))
    return nodes,outer,len(CHARTS)*(2650*nodes+650000068*outer+20)+1000
