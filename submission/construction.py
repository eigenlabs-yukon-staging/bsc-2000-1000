"""Code and RAM schedule for the aligned direct-address product."""
from collections import deque

MAIN_COLUMNS = (
    61033, 6582, 38677, 30442, 60646, 27952, 29717, 9118,
    3391, 33135, 21142, 52106, 11993, 56926, 29492, 54445,
)
TAIL_COLUMNS = (116, 202, 30, 61, 167, 141, 91, 242)
PLANE_GROUPS = 3
MAIN_BLOCKS = (992 - 256 * PLANE_GROUPS) // 16
TAIL_BASE = (1 << 32) - (1 << 16)
REPLICATE_LANES = sum(1 << i for i in (0, 16, 32, 48))


def parity_table(columns):
    values = [0]
    for col in columns:
        values.extend(v ^ col for v in values[:])
    return values


def leaders(columns):
    """Minimum-weight syndrome leaders, breaking ties by BFS discovery."""
    width = len(columns)
    syndrome_cols = columns + tuple(1 << i for i in range(width))
    size = 1 << width
    chosen = [-1] * size
    chosen[0] = 0
    queue = deque([0])
    while queue:
        s = queue.popleft()
        for i, col in enumerate(syndrome_cols):
            neighbor = s ^ col
            if chosen[neighbor] == -1:
                chosen[neighbor] = chosen[s] ^ (1 << i)
                queue.append(neighbor)
    assert all(v >= 0 for v in chosen)
    return chosen


def corrected_data(block, width, parity, chosen):
    low = (1 << width) - 1
    a = block & low
    s = parity[a] ^ (block >> width)
    return (a ^ chosen[s]) & low


def table_word(address, main_parity, main_leaders, tail_parity, tail_leaders):
    assert 48 <= address < (1 << 32)
    if address < TAIL_BASE:
        return REPLICATE_LANES * corrected_data(address, 16, main_parity, main_leaders)
    return corrected_data(address - TAIL_BASE, 8, tail_parity, tail_leaders) << 32


def encode(message):
    result = 0
    word_mask = (1 << 64) - 1
    for plane in range(PLANE_GROUPS):
        words = [(message >> (256 * plane + 64 * j)) & word_mask for j in range(4)]
        x0, x1, x2, d = words
        a, b, c = x0 ^ d, x1 ^ d, x2 ^ d
        words = [a, b, c, d, a ^ b ^ d, a ^ c ^ d, b ^ c ^ d]
        for j, word in enumerate(words):
            result |= word << (64 * (7 * plane + j))
    main_parity, tail_parity = map(parity_table, (MAIN_COLUMNS, TAIL_COLUMNS))
    for j in range(MAIN_BLOCKS):
        a = (message >> (256 * PLANE_GROUPS + 16 * j)) & 65535
        result |= (a | (main_parity[a] << 16)) << (448 * PLANE_GROUPS + 32 * j)
    a = message >> 992
    return result | ((a | (tail_parity[a] << 8)) << 1984)


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
    for plane in range(PLANE_GROUPS):
        yield from hamming_schedule(plane)
    for j in range(MAIN_BLOCKS // 2):
        yield "load", "Y", 7 * PLANE_GROUPS + j
        for high in (False, True):
            if high:
                yield "shr", "R", "Y", 32
            else:
                yield "and_", "R", "Y", (1 << 32) - 1
            yield "load", "D", "R"
            word, lane = divmod(2 * j + high, 4)
            w = 4 * PLANE_GROUPS + word
            yield "and_", "D", "D", 65535 << (16 * lane)
            if lane == 0:
                yield "move", f"U{w}", "D"
            else:
                yield "or_", f"U{w}", f"U{w}", "D"
    yield "load", "R", 31
    yield "or_", "R", "R", TAIL_BASE
    yield "load", "D", "R"
    yield "or_", "U15", "U15", "D"
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
