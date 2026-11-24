"""Field-cube code, constants, preloaded linear maps, and decoder sizes."""
from math import comb

WIDTH = 1000
RADIUS_FORWARD = 112
RADIUS_REVERSE = 109
MAX_RADIUS = max(RADIUS_FORWARD, RADIUS_REVERSE)
MOD_LOW = 1 | (1 << 34) | (1 << 45) | (1 << 53)
MODULUS = (1 << WIDTH) | MOD_LOW
MASK = (1 << WIDTH) - 1
THIRD = ((1 << WIDTH) - 1) // 3
EXPONENT = (THIRD + 1) // 3
CHAIN = (0, 1, 0, 0, 1, 1, 0)

ARRAYS = ("YS", "YP", "KEY", "SQ", "CUBE", "BEST", "Z2", "Z4", "T", "G", "H", "J",
          "ROOT", "ROTATED", "P1", "P2")
BASE = {name: 48 + 16 * j for j, name in enumerate(ARRAYS)}
NEXT = 48 + 16 * len(ARRAYS)
LAST = NEXT + MAX_RADIUS + 1
MULT = LAST + MAX_RADIUS + 1
RAW = MULT + 1024
HISTORY = RAW + 32
REDUCE_ONE = HISTORY + 32 * (MAX_RADIUS + 1)
REDUCE_SIX = REDUCE_ONE + 2
DIRECTION_TABLE = REDUCE_SIX + 64
MAP_BASE = DIRECTION_TABLE + 32 * WIDTH


def map_tags():
    shifts = {1, 2, 4}
    count = 1
    for bit in CHAIN:
        shifts.add(6 * count)
        if bit:
            shifts.add(12 * count)
        count = 2 * count + bit
    return tuple(sorted(shifts)) + ("zeta", "zeta2")


MEMORY_END = MAP_BASE + len(map_tags()) * 125 * 4096


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


def product(a, b):
    return remainder(carryless(a, b))


def square(a):
    return product(a, a)


def cube(a):
    return product(a, square(a))


def encode(a):
    return a | (cube(a) << WIDTH)


def frobenius(a, steps):
    for _ in range(steps):
        a = square(a)
    return a


def power(a, exponent):
    result = 1
    for bit in bin(exponent)[2:]:
        result = square(result)
        if bit == "1":
            result = product(result, a)
    return result


def root_of_unity():
    return power(7, THIRD)


def inverse_chain(z):
    z2 = square(z)
    z4 = frobenius(z, 2)
    t = frobenius(product(product(z2, z4), z), 4)
    g = t
    count = 1
    for bit in CHAIN:
        g = product(g, frobenius(g, 6 * count))
        if bit:
            g = product(g, frobenius(t, 12 * count))
        count = 2 * count + bit
    return product(g, z2)


def reduction_table(shift):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << shift))


def map_columns(tag):
    if isinstance(tag, int):
        base = frobenius(2, tag)
        value = 1
    else:
        base = 2
        value = root_of_unity()
        if tag == "zeta2":
            value = square(value)
    for _ in range(WIDTH):
        yield value
        value = product(value, base)


def map_words(tag):
    """Words in address order, for one linear map."""
    columns = iter(map_columns(tag))
    for _ in range(125):
        entries = [0]
        for direction in (next(columns) for _ in range(8)):
            entries.extend(direction ^ a for a in entries[:])
        for value in entries:
            for h in range(16):
                yield (value >> (64 * h)) & ((1 << 64) - 1)


def direction_words():
    for i in range(WIDTH):
        e = 1 << i
        for value in (square(e), cube(e)):
            for h in range(16):
                yield (value >> (64 * h)) & ((1 << 64) - 1)


def operation_bound():
    counts = tuple(
        sum(comb(WIDTH, i) for i in range(radius + 1))
        for radius in (RADIUS_FORWARD, RADIUS_REVERSE)
    )
    square_ops = 6641
    windows = (WIDTH + 5) // 6
    crossings = sum(6 * j % 64 > 58 for j in range(windows))
    multiply_ops = 16 + 32 + 5 * (51 + 16) + 57 * 64 + windows * (5 + 51 + 48) + 3 * crossings + 16
    forward = 293 + 21 + 98 + 3 * 562 + 226 + 13 + 98 + 6
    reverse = 14 * (square_ops + multiply_ops) + 32 + 66 + 3 * 293 + 2 * square_ops + 40
    return counts, forward * counts[0] + reverse * counts[1] + 227 + square_ops + multiply_ops
