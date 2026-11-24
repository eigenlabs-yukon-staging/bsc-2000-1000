"""Coupled field-cube encoder, fixed tables and the two-search operation bound."""
from math import comb

M = 500
WORDS = 8
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 8) | (1 << 13) | (1 << 21)
MODULUS = (1 << M) | MOD_LOW
BANDS = ((0, 60, 206), (61, 64, 194), (65, 68, 182), (69, 72, 170))
SLICE = 14
SLICES = (M + SLICE - 1) // SLICE
SLICE_ENTRIES = 1 << SLICE
AFFINE_STRIDE = SLICES * SLICE_ENTRIES * WORDS
RADIUS = BANDS[-1][1]
SKIP = M + 1
MASK64 = (1 << 64) - 1
SQUARE_BYTES = tuple(
    sum(((value >> bit) & 1) << (2 * bit) for bit in range(8)).to_bytes(2, "little")
    for value in range(256)
)

FIELD_ARRAYS = ("YA", "YB", "YC", "YD", "OUTER", "INNER", "OCUBE", "ICUBE", "BASECUBE", "SQ", "BEST_B", "BEST_C")
ARRAY = {name: 48 + WORDS * i for i, name in enumerate(FIELD_ARRAYS)}
TRAVERSAL_BASE = 48 + WORDS * len(FIELD_ARRAYS)
STACK_NEXT = tuple(TRAVERSAL_BASE + 2 * (M + 1) * i for i in range(3))
STACK_LAST = tuple(base + M + 1 for base in STACK_NEXT)
DISAGREE = TRAVERSAL_BASE + 6 * (M + 1)
AGREE = DISAGREE + M
SAVE_D = AGREE + M
SAVE_A = SAVE_D + (M + 1) * WORDS
MULT = SAVE_A + (M + 1) * WORDS
REDUCE_1 = MULT + 16 * WORDS
REDUCE_4 = REDUCE_1 + 2
SQUARE_BASE = REDUCE_4 + 16
DEPTH_BASE = SQUARE_BASE + BYTES * 256 * WORDS
CAP_BASE = DEPTH_BASE + RADIUS + 1
AFFINE_BASE = CAP_BASE + len(BANDS) * (M + 1)
MEMORY_END = AFFINE_BASE + M * AFFINE_STRIDE


def carryless(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def poly_remainder(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def remainder(value):
    while value >> M:
        upper = value >> M
        value = (value & MASK) ^ upper ^ (upper << 8) ^ (upper << 13) ^ (upper << 21)
    return value


def square(value):
    raw = value.to_bytes((value.bit_length() + 7) // 8, "little")
    return remainder(int.from_bytes(b"".join(SQUARE_BYTES[byte] for byte in raw), "little"))


def cube(value):
    return remainder(carryless(value, square(value)))


def codeword(b, c):
    return c ^ cube(b), b, c, b ^ cube(c)


def encode(message):
    return sum(value << (M * i) for i, value in enumerate(codeword(message & MASK, message >> M)))


def square_entry(byte, value):
    return square(value << (8 * byte))


def reduction_entry(c):
    return carryless(c, MOD_LOW)


def inner_cap(disagreements, budget):
    if disagreements > budget:
        return -1
    remaining = M - disagreements
    limit = 1 << (budget - disagreements)
    volume = 0
    for r in range(remaining + 1):
        volume += comb(remaining, r)
        if volume > limit:
            return r - 1
    return remaining


def affine_entry(direction, index, key):
    e = 1 << direction
    result = cube(e) if index == 0 else 0
    for bit in range(SLICE):
        j = index * SLICE + bit
        if j < M and (key >> bit) & 1:
            x = 1 << j
            result ^= remainder((square(e) << j) ^ (square(x) << direction))
    return result


def table_word(address):
    if REDUCE_1 <= address < REDUCE_4:
        return reduction_entry(address - REDUCE_1)
    if REDUCE_4 <= address < SQUARE_BASE:
        return reduction_entry(address - REDUCE_4)
    if SQUARE_BASE <= address < DEPTH_BASE:
        index, word = divmod(address - SQUARE_BASE, WORDS)
        byte, value = divmod(index, 256)
        return (square_entry(byte, value) >> (64 * word)) & MASK64
    if DEPTH_BASE <= address < CAP_BASE:
        depth = address - DEPTH_BASE
        band = next(i for i, (low, high, _) in enumerate(BANDS) if low <= depth <= high)
        return CAP_BASE + (M + 1) * band
    if CAP_BASE <= address < AFFINE_BASE:
        band, disagree = divmod(address - CAP_BASE, M + 1)
        radius = inner_cap(disagree, BANDS[band][2])
        return SKIP if radius < 0 else radius
    if AFFINE_BASE <= address < MEMORY_END:
        direction, offset = divmod(address - AFFINE_BASE, AFFINE_STRIDE)
        row, word = divmod(offset, WORDS)
        index, key = divmod(row, SLICE_ENTRIES)
        return (affine_entry(direction,index,key) >> (64*word)) & MASK64
    raise ValueError("not a fixed table address")


def operation_bound():
    k = WORDS
    shift = 3 * k + 3
    sq = BYTES * (5 + 3 * k) + k
    mul = 3 * k + 3 * (shift + k) + 44 * k + NIBBLES * (5 + shift + 3 * k) + k
    outer_nodes = sum(comb(M, i) for i in range(RADIUS + 1))
    inner_budget = sum(
        sum(comb(M, i) for i in range(low, high + 1)) << power
        for low, high, power in BANDS
    )
    outer_visit = sq + mul + 14 * M + 6 + 4 * k + 1
    crosses = sum((j*SLICE)%64+SLICE > 64 for j in range(SLICES))
    update = 2 + 2*k + SLICES*(6+3*k) + 3*crosses
    save_restore = 2*(2+3*k)
    inner_visit = 22*k + 2 + 44 + 45 + update + save_restore
    fixed = 100000 + 2 * (2*k + 1 + sq + mul)
    W = fixed + 2 * ((outer_visit + 44) * outer_nodes + inner_visit * inner_budget)
    return W, (sq, mul, outer_nodes, inner_budget, outer_visit, inner_visit, fixed, update, save_restore)
