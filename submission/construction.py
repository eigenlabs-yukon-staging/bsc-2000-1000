"""Explicit 64-bit components, address ranges, and the fixed RAM schedule."""
from collections import deque

GENERATOR = 141459731
EXTRA = (30, 27, 30, 27, 1, 2, 2, 11, 26, 5, 23, 25, 21, 27, 9, 8, 19, 6, 19, 1, 18, 21, 5, 31, 13, 20, 12, 25, 23, 27, 31, 16)
TAIL_COLUMNS = (116, 202, 30, 61, 167, 141, 91, 242)
OFFSET = 0x579a3c41
LOW_BASE = (1 << 32) - 3 * (1 << 16)
HIGH_BASE = LOW_BASE + (1 << 16)
TAIL_BASE = HIGH_BASE + (1 << 16)
MASK32 = (1 << 32) - 1
MASK16 = (1 << 16) - 1


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
    if 48 <= address < LOW_BASE:
        return chosen_main[address ^ OFFSET] & MASK32
    if LOW_BASE <= address < HIGH_BASE:
        return parity(address - LOW_BASE, MAIN_COLUMNS[:16]) ^ OFFSET
    if HIGH_BASE <= address < TAIL_BASE:
        return parity(address - HIGH_BASE, MAIN_COLUMNS[16:])
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
        yield "shr", "S", "Y", 32
        yield "and_", "L", "A", MASK16
        yield "shr", "H", "A", 16
        yield "or_", "L", "L", LOW_BASE
        yield "load", "L", "L"
        yield "or_", "H", "H", HIGH_BASE
        yield "load", "H", "H"
        yield "xor", "S", "S", "L"
        yield "xor", "S", "S", "H"
        yield "load", "C", "S"
        yield "xor", "A", "A", "C"
        dest = f"U{j//2}"
        if j % 2:
            yield "shl", "A", "A", 32
            yield "or_", dest, dest, "A"
        else:
            yield "move", dest, "A"
    yield "load", "Y", 31
    yield "or_", "Y", "Y", TAIL_BASE
    yield "load", "C", "Y"
    yield "or_", "U15", "U15", "C"
    yield "and_", "U15", "U15", (1 << 40) - 1
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
