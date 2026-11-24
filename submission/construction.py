"""Code and RAM schedule for the aligned direct-address product."""
from collections import deque

MAIN_COLUMNS = (
    61033, 6582, 38677, 30442, 60646, 27952, 9118, 3391, 33135, 52106,
    11993, 56926, 29492, 54445,
)
TAG_OFFSETS = (0,29717,25841,4324)
LOW_MASK = (1 << 14)-1
TAIL_COLUMNS = (116, 202, 30, 61, 167, 141, 91, 242)
TAIL_BASE = (1 << 32) - (1 << 16)


def parity_table(columns):
    values = [0]
    for col in columns:
        values.extend(v ^ col for v in values[:])
    return values


def leaders(columns, width=None):
    """Minimum-weight syndrome leaders, breaking ties by BFS discovery."""
    width = len(columns) if width is None else width
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


def quotient_choice(s, tag, leader):
    return min((leader[s^offset].bit_count()+(tag^t).bit_count(),t)
               for t,offset in enumerate(TAG_OFFSETS))


def main_decode(block, parity, leader):
    low=block & LOW_MASK
    tag=(block >> 14) & 3
    syndrome=(block >> 16)^parity[low]
    _, chosen=quotient_choice(syndrome,tag,leader)
    error=leader[syndrome^TAG_OFFSETS[chosen]]
    return ((low^error)&LOW_MASK)|(chosen<<14)


def table_word(address, main_parity, main_leaders, tail_parity, tail_leaders):
    assert 48 <= address < (1 << 32)
    if address < TAIL_BASE:
        return main_decode(address, main_parity, main_leaders)
    return corrected_data(address - TAIL_BASE, 8, tail_parity, tail_leaders) << 32


def encode(message):
    result = 0
    main_parity, tail_parity = map(parity_table, (MAIN_COLUMNS, TAIL_COLUMNS))
    for j in range(62):
        a = (message >> (16 * j)) & 65535
        parity=main_parity[a&LOW_MASK]^TAG_OFFSETS[a>>14]
        result |= (a | (parity << 16)) << (32 * j)
    a = message >> 992
    return result | ((a | (tail_parity[a] << 8)) << 1984)


def schedule():
    for j in range(31):
        yield "load", "Y", j
        for high in (False, True):
            if high:
                yield "shr", "R", "Y", 32
            else:
                yield "and_", "R", "Y", (1 << 32) - 1
            yield "load", "D", "R"
            w, lane = divmod(2 * j + high, 4)
            if lane:
                yield "shl", "D", "D", 16*lane
            if lane == 0:
                yield "move", f"U{w}", "D"
            else:
                yield "or_", f"U{w}", f"U{w}", "D"
    yield "load", "R", 31
    yield "or_", "R", "R", TAIL_BASE
    yield "load", "D", "R"
    yield "or_", "U15", "U15", "D"
    yield "and_", "U15", "U15", (1 << 40)-1
    for j in range(16):
        yield "store", None, 32 + j, f"U{j}"
