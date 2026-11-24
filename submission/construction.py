"""Explicit field component, preloaded data, and counted sizes."""
from math import comb

WIDTH = 500
RADIUS = 69
COMPONENTS = 1000 // WIDTH
WORDS = (WIDTH + 63) // 64
HIGH_BITS = WIDTH - 64 * (WORDS - 1)
BYTES = (WIDTH + 7) // 8
NIBBLES = (WIDTH + 3) // 4
MOD_LOW = 1 | sum(1 << i for i in (8, 13, 21))
MODULUS = (1 << WIDTH) | MOD_LOW

ARRAYS = ("YS", "YP", "KEY", "SQ", "CUBE", "BEST")
BASE = {name: 48 + WORDS * j for j, name in enumerate(ARRAYS)}
NEXT = 48 + WORDS * len(ARRAYS)
LAST = NEXT + RADIUS + 1
MULT = LAST + RADIUS + 1
REDUCE_ONE = MULT + 16 * WORDS
REDUCE_FOUR = REDUCE_ONE + 2
SQUARE_BASE = REDUCE_FOUR + 16
SAVED = SQUARE_BASE + BYTES * 256 * WORDS
AFFINE = SAVED + (RADIUS + 1) * WORDS
CHUNK_BITS = 14
CHUNKS = (WIDTH + CHUNK_BITS - 1) // CHUNK_BITS
CHUNK_ENTRIES = 1 << CHUNK_BITS
AFFINE_STRIDE = CHUNKS * CHUNK_ENTRIES * WORDS
MEMORY_END = AFFINE + WIDTH * AFFINE_STRIDE


def carryless(a, b):
    result = 0
    while b:
        lo = b & -b
        result ^= a << (lo.bit_length() - 1)
        b ^= lo
    return result


def remainder(a, b=MODULUS):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def square(a):
    return remainder(carryless(a, a))


def cube(a):
    return remainder(carryless(a, square(a)))


def encode(message):
    result = 0
    for i in range(COMPONENTS):
        payload = (message >> (i * WIDTH)) & ((1 << WIDTH) - 1)
        result |= (payload | (cube(payload) << WIDTH)) << (2 * i * WIDTH)
    return result


def reduction_table(shift):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << shift))


def square_table():
    for j in range(BYTES):
        values = [0]
        for k in range(8):
            direction = square(1 << (8 * j + k))
            values.extend(direction ^ a for a in values[:])
        for value in values:
            for h in range(WORDS):
                yield (value >> (64 * h)) & ((1 << 64) - 1)


def derivative_columns(i):
    e=1 << i
    linear,quadratic=square(e),e
    for _ in range(WIDTH):
        yield linear ^ quadratic
        linear=remainder(linear << 1)
        quadratic=remainder(quadratic << 2)
    yield from (0 for _ in range(CHUNKS*CHUNK_BITS-WIDTH))


def affine_entry(i,chunk,key):
    columns=tuple(derivative_columns(i))
    result=cube(1 << i) if chunk == 0 else 0
    for j in range(CHUNK_BITS):
        if key >> j & 1:
            result ^= columns[CHUNK_BITS*chunk+j]
    return result


def operation_bound():
    k = WORDS
    shift = 3 * k + 3
    square_ops = BYTES * (5 + 3 * k) + k
    multiply_ops = 3 * k + 3 * (shift + k) + 44 * k + NIBBLES * (5 + shift + 3 * k) + k
    scoring_ops = 18 * k + 5
    nodes = sum(comb(WIDTH, w) for w in range(RADIUS + 1))
    crossings=sum((j*CHUNK_BITS)%64+CHUNK_BITS>64 for j in range(CHUNKS))
    update_ops=2+2*k+CHUNKS*(6+3*k)+3*crossings
    save_restore=2*(2+3*k)
    ops=COMPONENTS*((update_ops+save_restore+scoring_ops+40)*nodes
                    +24*k+1+square_ops+multiply_ops)
    return nodes, ops, (square_ops, multiply_ops, scoring_ops)
