"""Explicit hybrid components, shared table, and fixed RAM schedule."""
from itertools import combinations

FIRST_ROWS = ('0000001100111111', '0001011010101111')
CHANGED_PARITY_BIT = 27
CHANGED_DATA_MASK = 0x2e27aa03
OFFSET = 0x1739e3b1
TAIL_LIMIT = 1 << 24
MASK32 = (1 << 32) - 1


def row_bits(row):
    return sum(int(bit) << i for i, bit in enumerate(row))


def rotate16(value, i):
    return ((value << i) | (value >> (16 - i))) & 0xffff


def transpose_row(value):
    return sum(((value >> (-i % 16)) & 1) << i for i in range(16))


def main_columns():
    a, b = map(row_bits, FIRST_ROWS)
    rows = ((a, b), (transpose_row(b), transpose_row(a)))
    columns = [rotate16(lo, i) | (rotate16(hi, i) << 16)
               for lo, hi in rows for i in range(16)]
    return tuple(c ^ (((CHANGED_DATA_MASK >> i) & 1) << CHANGED_PARITY_BIT)
                 for i, c in enumerate(columns))


PARITY_ORDER = (14, 15, 25, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26)


def permute_parity(value):
    return sum(((value >> i) & 1) << j for j, i in enumerate(PARITY_ORDER))


BASE_COLUMNS = main_columns()
BASE_SYNDROME_COLUMNS = BASE_COLUMNS + tuple(1 << i for i in range(32))
MAIN_COLUMNS = tuple(permute_parity(c) for c in BASE_COLUMNS)
SYNDROME_COLUMNS = MAIN_COLUMNS + tuple(1 << i for i in range(32))


def parity(value, columns):
    result = 0
    for i, column in enumerate(columns):
        if value >> i & 1:
            result ^= column
    return result


TAIL_COLUMNS = tuple(c >> 8 for c in MAIN_COLUMNS)
TAIL_SYNDROME_COLUMNS = TAIL_COLUMNS + tuple(1 << i for i in range(24))
MASK40 = (1 << 40) - 1
MASK64 = (1 << 64) - 1


def encode(message):
    result = 0
    words = [(message >> (64 * j)) & MASK64 for j in range(4)]
    x0, x1, x2, d = words
    a, b, c = x0 ^ d, x1 ^ d, x2 ^ d
    words = [a, b, c, d, a ^ b ^ d, a ^ c ^ d, b ^ c ^ d]
    for j, word in enumerate(words):
        result |= word << (64 * j)
    for j in range(22):
        a = (message >> (256 + 32 * j)) & MASK32
        result |= (a | (parity(a, MAIN_COLUMNS) << 32)) << (64 * (7 + j))
    tail = message >> 960
    return result | ((tail | (parity(tail & MASK32, TAIL_COLUMNS) << 40)) << (64 * 29))


def minimum_leader(syndrome, columns):
    """Offline table specification, ordered first by weight then coordinates."""
    for weight in range(len(columns) + 1):
        for indices in combinations(range(len(columns)), weight):
            value = 0
            for i in indices:
                value ^= columns[i]
            if value == syndrome:
                return sum(1 << i for i in indices)
    raise ValueError("syndrome not spanned")


def table_word(address):
    """Emit one fixed word of the preloaded table."""
    if not 48 <= address < 1 << 32:
        raise ValueError("outside fixed table")
    check = parity(address, MAIN_COLUMNS) ^ OFFSET
    if address < TAIL_LIMIT:
        check &= ~255
        correction = minimum_leader(address ^ (OFFSET >> 8), TAIL_SYNDROME_COLUMNS) & MASK32
    else:
        correction = minimum_leader(address ^ OFFSET, SYNDROME_COLUMNS) & MASK32
    return (check << 32) | correction


def hamming_schedule(plane):
    for j, name in enumerate(("A", "B", "C", "D", "P0", "P1", "P2")):
        yield "load", name, 7 * plane + j
    yield "xor", "L0", "A", "D"
    yield "xor", "L1", "B", "D"
    yield "xor", "L2", "C", "D"
    yield "xor", "S0", "L0", "B"
    yield "xor", "S0", "S0", "P0"
    yield "xor", "S1", "L0", "C"
    yield "xor", "S1", "S1", "P1"
    yield "xor", "S2", "L1", "C"
    yield "xor", "S2", "S2", "P2"
    yield "and_", "M0", "S0", "S1"
    yield "and_", "M1", "S0", "S2"
    yield "and_", "M2", "S1", "S2"
    yield "and_", "M3", "M0", "S2"
    for j, name in enumerate(("L0", "L1", "L2", "D")):
        yield "xor", f"U{4 * plane + j}", name, f"M{j}"


def schedule():
    yield from hamming_schedule(0)
    for j in range(22):
        yield "load", "Y", 7 + j
        yield "and_", "A", "Y", MASK32
        yield "load", "Check", "A"
        yield "xor", "Marked", "Y", "Check"
        yield "shr", "S", "Marked", 32
        yield "load", "Correction", "S"
        yield "xor", "Decoded", "Y", "Correction"
        dest = f"U{4 + j // 2}"
        if j % 2:
            yield "shl", "Decoded", "Decoded", 32
            yield "or_", dest, dest, "Decoded"
        else:
            yield "and_", dest, "Decoded", MASK32
    yield "load", "Y", 29
    yield "and_", "A", "Y", MASK32
    yield "load", "Check", "A"
    yield "xor", "Marked", "Y", "Check"
    yield "shr", "S", "Marked", 40
    yield "load", "Correction", "S"
    yield "xor", "Decoded", "Y", "Correction"
    yield "and_", "U15", "Decoded", MASK40
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
