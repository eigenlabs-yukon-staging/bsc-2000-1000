"""The explicit four-group field code and its fixed tables."""
from functools import lru_cache
from itertools import combinations
from math import comb

M = 500
MASK = (1 << M) - 1
HALF = 250
SUBMASK = (1 << HALF) - 1
MOD_LOW = 1 | (1 << 32) | (1 << 40) | (1 << 59)
MODULUS = (1 << HALF) | MOD_LOW
XI = 1 << 191
THETA = 1 << HALF
LAMBDA = 2
THIRD = MASK // 3
CUBE_EXP = (2 * THIRD + 1) // 3
RADII = (105, 105, 105, 105, 105, 105)
SPLIT_RADIUS = 106
MAX_RADIUS = max(*RADII, SPLIT_RADIUS)
BYTES = 64
PAIRS = tuple(combinations(range(4), 2))
CHART_COSTS = (20000, 20000, 152000, 20000, 152000, 310000)
SPLIT_COST = 3000
SPLIT_CHARTS = tuple((outer, j, k) for outer in (0, 1)
                     for j in (0, 1) for k in (0, 1))
ARRAY_LIST = (
    "Y0", "Y1", "Y2", "Y3", "KEY0", "KEY1", "V0", "V1", "V2", "V3",
    "BEST0", "BEST1", "X2", "X4", "T", "G", "H", "J", "CHECK", "C1", "C2",
    "INV", "INV2", "TEMP", "K", "TARGET", "ASROOT", "A0", "A1", "BASEA", "BASEB",
    "ROOT", "INVA", "KNOWN", "M0", "M1", "M2", "M3", "M4",
    "RATIO", "INVCOEF", "DELTA", "F1", "F2", "F3", "LIN", "QUAD",

)
BASE = {name: 48 + 8 * i for i, name in enumerate(ARRAY_LIST)}
NEXT = 48 + 8 * len(ARRAY_LIST)
LAST = NEXT + MAX_RADIUS + 1
INEXT = LAST + MAX_RADIUS + 1
ILAST = INEXT + MAX_RADIUS + 1
MULT = ILAST + MAX_RADIUS + 1
REDSHIFT1 = MULT + 16 * 4
REDSHIFT4 = REDSHIFT1 + 2
XI_MAP = REDSHIFT4 + 16
MAP_BASE = XI_MAP + 32 * 256 * 4


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


def geom(a, length, step):
    result = a
    k = 1
    for bit in bin(length)[3:]:
        result = product(result, frobenius(result, step * k))
        if bit == "1":
            result = product(result, frobenius(a, 2 * step * k))
        k = 2 * k + int(bit)
    return result


def inverse(a):
    assert a
    return square(geom(a, M - 1, 1))


def root_power(a):
    a7 = product(product(square(a), frobenius(a, 2)), a)
    t = frobenius(a7, 3)
    v = product(a, geom(t, 83, 6))
    return v


def cube_root(a):
    v = root_power(a)
    return v if cube(v) == a else None


@lru_cache(None)
def constant(tag):
    if tag == "lambda":
        return LAMBDA
    if tag == "lambda_inv":
        return inverse(LAMBDA)
    if tag == "lambda_minus_e":
        return power(LAMBDA, MASK - CUBE_EXP)
    if tag == "zeta":
        return power(LAMBDA, THIRD)
    if tag == "zeta2":
        return square(constant("zeta"))
    if tag == "theta":
        return THETA
    if tag == "theta2":
        return square(THETA)
    raise ValueError("unknown field constant")


@lru_cache(None)
def AS_pivots():
    pivots = {}
    for i in range(M):
        preimage = 1 << i
        image = square(preimage) ^ preimage
        while image:
            pivot = image.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (image, preimage)
                break
            r, x = pivots[pivot]
            image ^= r
            preimage ^= x
    return tuple((i, *pivots[i]) for i in sorted(pivots, reverse=True))


def AS_map(y):
    x = 0
    for pivot, image, preimage in AS_pivots():
        if (y >> pivot) & 1:
            y ^= image
            x ^= preimage
    return x


def cube_roots(x):
    base = cube_root(x)
    if base is None:
        return ()
    return (base, product(constant("zeta"), base), product(constant("zeta2"), base))


def codeword(a, b):
    return a, b, product(a, b), cube(a) ^ product(LAMBDA, cube(b))


def encode(message):
    coordinates = codeword(message & MASK, message >> M)
    return sum(value << (M * i) for i, value in enumerate(coordinates))


def chart(pair, x, y):
    """Chart candidates passing the good-message test."""
    if pair == (2, 3) and not y:
        return ()
    if not x:
        return ()
    return tuple((a,b) for a,b in raw_chart(pair,x,y) if good(a) and good(b))


def raw_chart(pair, x, y):
    """Host semantic counterpart of the bounded field chart macros."""
    if pair == (0, 1):
        return ((x, y),) if x and y else ()
    if not x:
        return ()
    if pair == (0, 2):
        return ((x, product(y, inverse(x))),) if y else ()
    if pair == (1, 2):
        return ((product(y, inverse(x)), x),) if y else ()
    if pair == (0, 3):
        target = product(y ^ cube(x), constant("lambda_inv"))
        return tuple((x, b) for b in cube_roots(target)) if target else ()
    if pair == (1, 3):
        target = y ^ product(LAMBDA, cube(x))
        return tuple((a, x) for a in cube_roots(target)) if target else ()
    assert pair == (2, 3)
    c, d = x, y
    if not d:
        return ()
    k = product(LAMBDA, cube(c))
    z = product(k, square(inverse(d)))
    w = AS_map(z)
    if square(w) ^ w != z:
        return ()
    aa = product(d, w)
    root = root_power(aa)
    test = cube(root)
    if test == aa:
        a, b = root, product(c, inverse(root))
    elif test == product(constant("zeta2"), aa):
        b = product(root, constant("lambda_minus_e"))
        a = product(c, inverse(b))
    else:
        return ()
    return ((a, b),
            (product(constant("zeta"), a), product(constant("zeta2"), b)),
            (product(constant("zeta2"), a), product(constant("zeta"), b)))


def map_tags():
    shifts = {1, 2, 3}
    for length, step in ((83, 6), (M - 1, 1)):
        k = 1
        for bit in bin(length)[3:]:
            shifts.add(step * k)
            if bit == "1":
                shifts.add(2 * step * k)
            k = 2 * k + int(bit)
    return tuple(sorted(shifts)) + ("lambda", "lambda_inv", "lambda_minus_e", "zeta", "zeta2", "AS", "theta", "theta2")


def map_columns(tag):
    for half in (0, 1):
        for i in range(8 * (BYTES // 2)):
            if i >= HALF:
                yield 0
                continue
            e = 1 << (half * HALF + i)
            if tag == "AS":
                yield AS_map(e)
            elif isinstance(tag, int):
                yield frobenius(e, tag)
            else:
                yield product(constant(tag), e)


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


DIRECTIONS = MAP_BASE + len(map_tags()) * BYTES * 256 * 8
SAVED = DIRECTIONS + M * 16
AFFINE = SAVED + (MAX_RADIUS+1) * 24
AFFINE_STRIDE = BYTES * 256 * 8
MEMORY_END = AFFINE + M * AFFINE_STRIDE


def W_bound():
    radii = RADII + (SPLIT_RADIUS,) * len(SPLIT_CHARTS)
    costs = CHART_COSTS + (SPLIT_COST,) * len(SPLIT_CHARTS)
    nodes = tuple(sum(comb(2 * M, i) for i in range(r + 1)) for r in radii)
    outer = tuple(sum(comb(M, i) for i in range(r + 1)) for r in radii)
    w = sum((cost + 50) * n + 650000068 * o
            for cost, n, o in zip(costs, nodes, outer)) + 1000
    return nodes, outer, w
