"""Explicit hybrid components, shared table, and fixed RAM schedule."""
from itertools import combinations

GENERATOR = 141459731
EXTRA = (29, 24, 30, 9, 28, 26, 2, 0, 18, 9, 13, 12, 16, 11, 4, 23, 6, 8, 14, 0, 7, 19, 14, 28, 0, 3, 2, 31, 22, 15, 13, 0)
OFFSET = 0xf7c90230
TAIL_LIMIT = 1 << 24
MASK32 = (1 << 32) - 1


def remainder(value):
    while value.bit_length() > 27:
        value ^= GENERATOR << (value.bit_length() - 28)
    return value


def parity_column(i):
    shifted = (1 << (27 + i)) | (EXTRA[i] << 59)
    low_parity = remainder(shifted)
    return low_parity | (EXTRA[i] << 27)


MAIN_COLUMNS = tuple(parity_column(i) for i in range(32))
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
    for plane in range(3):
        words = [(message >> (256 * plane + 64 * j)) & MASK64 for j in range(4)]
        x0, x1, x2, d = words
        a, b, c = x0 ^ d, x1 ^ d, x2 ^ d
        words = [a, b, c, d, a ^ b ^ d, a ^ c ^ d, b ^ c ^ d]
        for j, word in enumerate(words):
            result |= word << (64 * (7 * plane + j))
    for j in range(6):
        a = (message >> (768 + 32 * j)) & MASK32
        result |= (a | (parity(a, MAIN_COLUMNS) << 32)) << (64 * (21 + j))
    tail = message >> 960
    return result | ((tail | (parity(tail & MASK32, TAIL_COLUMNS) << 40)) << (64 * 27))


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
    for plane in range(3):
        yield from hamming_schedule(plane)
    for j in range(6):
        yield "load", "Y", 21 + j
        yield "and_", "A", "Y", MASK32
        yield "load", "Check", "A"
        yield "xor", "Marked", "Y", "Check"
        yield "shr", "S", "Marked", 32
        yield "load", "Correction", "S"
        yield "xor", "Decoded", "Y", "Correction"
        dest = f"U{12 + j // 2}"
        if j % 2:
            yield "shl", "Decoded", "Decoded", 32
            yield "or_", dest, dest, "Decoded"
        else:
            yield "and_", dest, "Decoded", MASK32
    yield "load", "Y", 27
    yield "and_", "A", "Y", MASK32
    yield "load", "Check", "A"
    yield "xor", "Marked", "Y", "Check"
    yield "shr", "S", "Marked", 40
    yield "load", "Correction", "S"
    yield "xor", "Decoded", "Y", "Correction"
    yield "and_", "U15", "Decoded", MASK40
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
