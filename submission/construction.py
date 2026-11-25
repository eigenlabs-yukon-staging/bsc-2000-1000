"""Explicit 64-bit components, address ranges, and the fixed RAM schedule."""
from collections import deque

FIRST_ROWS = ('0000001100111111', '0001011010101111')
CHANGED_PARITY_BIT = 27
CHANGED_DATA_MASK = 0x2e27aa03
TAIL_COLUMNS = (116, 202, 30, 61, 167, 141, 91, 242)
OFFSET = 0xa2aa23b1
TAIL_BASE = (1 << 32) - (1 << 16)
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


MAIN_COLUMNS = main_columns()
SYNDROME_COLUMNS = MAIN_COLUMNS + tuple(1 << i for i in range(32))


def parity(value, columns):
    result = 0
    for i, column in enumerate(columns):
        if value >> i & 1:
            result ^= column
    return result


def tail_leaders():
    chosen = [-1] * 256
    chosen[0] = 0
    queue = deque([0])
    columns = TAIL_COLUMNS + tuple(1 << i for i in range(8))
    while queue:
        syndrome = queue.popleft()
        for i, col in enumerate(columns):
            nxt = syndrome ^ col
            if chosen[nxt] == -1:
                chosen[nxt] = chosen[syndrome] ^ (1 << i)
                queue.append(nxt)
    return tuple(chosen)


def encode(message):
    result = 0
    for i in range(31):
        a = (message >> (32 * i)) & MASK32
        result |= (a | (parity(a, MAIN_COLUMNS) << 32)) << (64 * i)
    a = message >> 992
    return result | ((a | (parity(a, TAIL_COLUMNS) << 8)) << 1984)


def table_word(address, chosen_main, chosen_tail):
    if 48 <= address < TAIL_BASE:
        check = parity(address, MAIN_COLUMNS) ^ OFFSET
        correction = chosen_main[address ^ OFFSET] & MASK32
        return (check << 32) | correction
    if TAIL_BASE <= address < 1 << 32:
        r = address - TAIL_BASE
        a = r & 255
        syndrome = (r >> 8) ^ parity(a, TAIL_COLUMNS)
        return ((a ^ chosen_tail[syndrome]) & 255) << 32
    raise ValueError("outside fixed tables")


def schedule():
    for j in range(31):
        yield "load", "Y", j
        yield "and_", "A", "Y", MASK32
        yield "load", "Check", "A"
        yield "xor", "Marked", "Y", "Check"
        yield "shr", "S", "Marked", 32
        yield "load", "Correction", "S"
        yield "xor", "Decoded", "Y", "Correction"
        dest = f"U{j//2}"
        if j % 2:
            yield "shl", "Decoded", "Decoded", 32
            yield "or_", dest, dest, "Decoded"
        else:
            yield "and_", dest, "Decoded", MASK32
    yield "load", "Y", 31
    yield "or_", "Y", "Y", TAIL_BASE
    yield "load", "C", "Y"
    yield "or_", "U15", "U15", "C"
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
