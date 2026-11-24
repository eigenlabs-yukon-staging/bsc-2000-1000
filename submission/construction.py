"""The explicit field code and decoder constants."""
from math import comb

WIDTH = 1000
RADIUS = 130
MOD_LOW = 1 | (1 << 34) | (1 << 45) | (1 << 53)
MODULUS = (1 << WIDTH) | MOD_LOW
MASK = (1 << WIDTH) - 1

YS = 48
YP = YS + 16
U = YP + 16
SQ = U + 16
CUBE = SQ + 16
BEST = CUBE + 16
NEXT = BEST + 16
LAST = NEXT + RADIUS + 1
MULT = LAST + RADIUS + 1
SQUARE_TABLE = MULT + 256
REDUCE_ONE = SQUARE_TABLE + 125 * 4096
REDUCE_FOUR = REDUCE_ONE + 2
MEMORY_END = REDUCE_FOUR + 16


def remainder(a, b=MODULUS):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def carryless(a, b):
    result = 0
    while b:
        lo = b & -b
        result ^= a << (lo.bit_length() - 1)
        b ^= lo
    return result


def square(a):
    return remainder(carryless(a, a))


def cube(a):
    return remainder(carryless(a, square(a)))


def encode(a):
    return a | (cube(a) << WIDTH)


def square_table():
    """Words in address order, one 16-word field element for each byte."""
    mask_word = (1 << 64) - 1
    for j in range(125):
        basis = [square(1 << (8 * j + k)) for k in range(8)]
        elements = [0]
        for direction in basis:
            elements.extend(direction ^ a for a in elements[:])
        for value in elements:
            for k in range(16):
                yield (value >> (64 * k)) & mask_word


def reduction_table(shift):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << shift))


def operation_bound():
    candidates = sum(comb(WIDTH, r) for r in range(RADIUS + 1))
    return candidates, 33943 * candidates + 178
