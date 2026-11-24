"""Field definitions, table data and counted dimensions for the three-component code."""
from dataclasses import dataclass
from math import comb

ARRAYS = ("YS", "YP", "KEY", "SQ", "CUBE", "BEST")
CHUNK_BITS = 15
CHUNK_ENTRIES = 1 << CHUNK_BITS
RADIUS = 47


def carryless(a, b):
    result = 0
    while b:
        lo = b & -b
        result ^= a << (lo.bit_length() - 1)
        b ^= lo
    return result


def polynomial_remainder(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


@dataclass
class Field:
    width: int
    exponents: tuple
    rabin_degrees: tuple
    start: int

    def __post_init__(self):
        self.words = (self.width + 63) // 64
        self.high_bits = self.width - 64 * (self.words - 1)
        self.bytes = (self.width + 7) // 8
        self.nibbles = (self.width + 3) // 4
        self.chunks = (self.width + CHUNK_BITS - 1) // CHUNK_BITS
        self.low = 1 | sum(1 << i for i in self.exponents)
        self.modulus = (1 << self.width) | self.low
        self.base = {name: self.start + self.words * j for j, name in enumerate(ARRAYS)}
        self.next = self.start + self.words * len(ARRAYS)
        self.last = self.next + RADIUS + 1
        self.mult = self.last + RADIUS + 1
        self.reduce_one = self.mult + 16 * self.words
        self.reduce_four = self.reduce_one + 2
        self.square_base = self.reduce_four + 16
        self.saved = self.square_base + self.bytes * 256 * self.words
        self.affine = self.saved + (RADIUS + 1) * self.words
        self.affine_stride = self.chunks * CHUNK_ENTRIES * self.words
        self.end = self.affine + self.width * self.affine_stride

    def remainder(self, value):
        return polynomial_remainder(value, self.modulus)

    def square(self, a):
        return self.remainder(carryless(a, a))

    def cube(self, a):
        return self.remainder(carryless(a, self.square(a)))

    def reduction_table(self, shift):
        return tuple(carryless(c, self.low) for c in range(1 << shift))

    def square_entry(self, position, key):
        return self.square(key << (8 * position))

    def derivative_columns(self, direction):
        e = 1 << direction
        linear, quadratic = self.square(e), e
        for _ in range(self.width):
            yield linear ^ quadratic
            linear = self.remainder(linear << 1)
            quadratic = self.remainder(quadratic << 2)
        yield from (0 for _ in range(self.chunks * CHUNK_BITS - self.width))

    def affine_entry(self, direction, chunk, key):
        columns = tuple(self.derivative_columns(direction))
        result = self.cube(1 << direction) if chunk == 0 else 0
        for j in range(CHUNK_BITS):
            if key >> j & 1:
                result ^= columns[CHUNK_BITS * chunk + j]
        return result

    def operation_bound(self):
        k = self.words
        shift = 3 * k + 3
        square_ops = self.bytes * (5 + 3 * k) + k
        multiply_ops = 3 * k + 3 * (shift + k) + 44 * k + self.nibbles * (5 + shift + 3 * k) + k
        scoring_ops = 18 * k + 5
        nodes = sum(comb(self.width, weight) for weight in range(RADIUS + 1))
        crossings = sum((j * CHUNK_BITS) % 64 + CHUNK_BITS > 64 for j in range(self.chunks))
        update_ops = 2 + 2 * k + self.chunks * (6 + 3 * k) + 3 * crossings
        save_restore = 2 * (2 + 3 * k)
        ops = ((update_ops + save_restore + scoring_ops + 40) * nodes
               + 24 * k + 1 + square_ops + multiply_ops)
        return ops, (nodes, square_ops, multiply_ops, scoring_ops, update_ops, save_restore)


F333 = Field(333, (12, 47, 50), (9, 111), 48)
F334 = Field(334, (2, 9, 47), (2, 167), F333.end)
FIELDS = (F333, F333, F334)


def encode(message):
    codeword = 0
    position = 0
    for field in FIELDS:
        payload = (message >> position) & ((1 << field.width) - 1)
        block = payload | (field.cube(payload) << field.width)
        codeword |= block << (2 * position)
        position += field.width
    return codeword


def operation_bound():
    return sum(field.operation_bound()[0] for field in FIELDS)
