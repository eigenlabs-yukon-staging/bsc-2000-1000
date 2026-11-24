"""Field definitions, table data and counted dimensions for the bidirectional field-cube product."""
from dataclasses import dataclass
from math import comb

SQUARE_BYTES = tuple(
    sum(((value >> bit) & 1) << (2 * bit) for bit in range(8)).to_bytes(2, "little")
    for value in range(256)
)

ARRAYS = ("YS", "YP", "KEY", "SQ", "CUBE", "BEST", "T0", "T1", "T2", "Y", "Z", "CHECK")
CHUNK_BITS = 14
CHUNK_ENTRIES = 1 << CHUNK_BITS
RADIUS = 56
PARITY_RADIUS = 54
WIDTH = 500
EVEN = WIDTH % 2 == 0
DEGREE_MOD_SIX = WIDTH % 6
ROOT_COUNT = 3 if EVEN else 1
MAP_BITS = 16
MAP_ENTRIES = 1 << MAP_BITS


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


def root_program():
    """Fixed macro circuit from the parity candidate KEY to a cube root Z."""
    if EVEN:
        if DEGREE_MOD_SIX == 4:
            seed = (4, 5, 6)
            suffix_square = True
        else:
            assert DEGREE_MOD_SIX == 2
            seed = (3, 4, 5)
            suffix_square = False
        for i, exponent in enumerate(seed):
            yield ("MAP", exponent, "KEY", f"T{i}")
        yield ("MUL", "T0", "T1", "Y")
        yield ("MUL", "Y", "T2", "Y")
        base = "Y"
        stride = 6
        count = WIDTH // 6
    else:
        base = "KEY"
        stride = 2
        count = (WIDTH + 1) // 2
    yield ("COPY", base, "Z")
    j = 1
    for bit in bin(count)[3:]:
        yield ("MAP", stride * j, "Z", "T0")
        yield ("MUL", "Z", "T0", "Z")
        j *= 2
        if bit == "1":
            yield ("MAP", stride * j, base, "T0")
            yield ("MUL", "Z", "T0", "Z")
            j += 1
    if EVEN:
        if suffix_square:
            yield ("MAP", 1, "KEY", "T0")
            yield ("MUL", "Z", "T0", "Z")
        else:
            yield ("MUL", "Z", "KEY", "Z")


ROOT_PROGRAM = tuple(root_program())
MAPS = tuple(sorted({term[1] for term in ROOT_PROGRAM if term[0] == "MAP"} | ({1} if EVEN else set())))
if EVEN:
    MAPS += ("omega",)


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
        self.map_base = self.affine + self.width * self.affine_stride
        self.map_slices = (self.width + MAP_BITS - 1) // MAP_BITS
        self.map_stride = self.map_slices * MAP_ENTRIES * self.words
        self.end = self.map_base + len(MAPS) * self.map_stride

    def remainder(self, value):
        mask = (1 << self.width) - 1
        while value >> self.width:
            high = value >> self.width
            value = (value & mask) ^ high
            for i in self.exponents:
                value ^= high << i
        return value

    def square(self, a):
        octets = a.to_bytes((a.bit_length() + 7) // 8, "little")
        even = b"".join(SQUARE_BYTES[byte] for byte in octets)
        return self.remainder(int.from_bytes(even, "little"))

    def cube(self, a):
        return self.remainder(carryless(a, self.square(a)))

    def power(self, a, exponent):
        result = 1
        for bit in bin(exponent)[2:]:
            result = self.square(result)
            if bit == "1":
                result = self.remainder(carryless(result, a))
        return result

    def map_entry(self, tag, position, key):
        value = key << (MAP_BITS * position)
        if tag == "omega":
            return self.remainder(carryless(value, OMEGA))
        for _ in range(tag):
            value = self.square(value)
        return value

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
        map_ops = self.map_slices * (5 + 3 * k) + k
        root_maps = sum(instruction[0] == "MAP" for instruction in ROOT_PROGRAM)
        root_muls = sum(instruction[0] == "MUL" for instruction in ROOT_PROGRAM)
        root_ops = 2 * k + root_maps * map_ops + root_muls * multiply_ops
        if EVEN:
            root_ops += 2 * map_ops + multiply_ops + 8 * k + 2
        parity_nodes = sum(comb(self.width, w) for w in range(PARITY_RADIUS + 1))
        ops += 2 * k + 1 + (root_ops + ROOT_COUNT * scoring_ops + 40) * parity_nodes
        return ops, (nodes, parity_nodes, square_ops, multiply_ops, scoring_ops,
                     update_ops, save_restore, map_ops, root_ops)


FIELD = Field(500, (8, 13, 21), (100, 250), 48)
FIELDS = (FIELD,) * 2
OMEGA = 2807832145868348318082146064684578644799627959147856813580105169703142315075359039909955523629982985396826262832967958662726040844777023019727190692741


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
