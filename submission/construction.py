"""Fixed linear code, table functions and primitive count."""
from collections import deque

PRIVATE_BITS = 19
TAG_BITS = 12
PARITY_BITS = 31
SIDE_BITS = 50
BLOCKS = 10
GROUPS = 1 << (TAG_BITS - 3)
SLOTS = 1 << TAG_BITS
LEADER_BASE = 1 << 31
PARITY_BASE = 1 << 30
TAG_BASE = PARITY_BASE + (1 << PRIVATE_BITS)
GAMMA_BASE = TAG_BASE + SLOTS
PROFILE_BASE = GAMMA_BASE + SLOTS
PACKED_BASE = PROFILE_BASE + 4 * SLOTS
END = PACKED_BASE + (1 + 3 * 8) * GROUPS
HIGH = 0x8080808080808080
REPLICATE = 0x0101010101010101

# Images of the 19 private basis bits, followed by the 12 tag basis bits.
COLUMNS = (1051358863, 590810211, 1475472909, 1978738269, 651388888, 232409661, 1369593774, 2076077886, 429576749, 416644663, 1299312799, 1275935853, 711053259, 784118576, 130579930, 548027298, 1864563713, 618131721, 451792847, 344514561, 689029122, 705654788, 1815220232, 541642768, 2025672736, 1230393408, 1387044992, 1700348160, 849244672, 1561334784, 975185920)


def linear_map(value, columns):
    result = 0
    while value:
        bit = value & -value
        result ^= columns[bit.bit_length() - 1]
        value ^= bit
    return result


def parity(private):
    return linear_map(private, COLUMNS[:PRIVATE_BITS])


def displacement(tag):
    return linear_map(tag, COLUMNS[PRIVATE_BITS:])


def gamma(tag):
    even = (SLOTS - 1) // 3
    low = tag & even
    high = (tag >> 1) & even
    return high | ((low ^ high) << 1)


def tag_tuple(a, s):
    d = gamma(s)
    return a, a ^ s, a ^ d, a ^ s ^ d


def encode(message):
    word = 0
    private_mask = (1 << PRIVATE_BITS) - 1
    for block in range(BLOCKS):
        payload = (message >> (100 * block)) & ((1 << 100) - 1)
        tags = tag_tuple((payload >> 76) & (SLOTS - 1), payload >> 88)
        for side, tag in enumerate(tags):
            u = (payload >> (PRIVATE_BITS * side)) & private_mask
            word |= (u | ((parity(u) ^ displacement(tag)) << PRIVATE_BITS)) << (200 * block + SIDE_BITS * side)
    return word


def profile_address(side, tag):
    return PROFILE_BASE + SLOTS * side + tag


def packed_address(side, permutation, group):
    start = 0 if side == 0 else 1 + 8 * (side - 1)
    return PACKED_BASE + (start + permutation) * GROUPS + group


def directions():
    for s in range(1, SLOTS):
        d = gamma(s)
        if s < min(d, s ^ d):
            yield s, d


def syndrome_columns():
    return COLUMNS[:PRIVATE_BITS] + tuple(1 << i for i in range(PARITY_BITS))


def leaders(columns=None):
    if columns is None:
        columns = syndrome_columns()
    chosen = {0: 0}
    pending = deque([0])
    while pending:
        value = pending.popleft()
        for i, column in enumerate(columns):
            successor = value ^ column
            if successor not in chosen:
                chosen[successor] = chosen[value] ^ (1 << i)
                pending.append(successor)
    return chosen


def table_word(address, chosen):
    if LEADER_BASE <= address < 1 << 32:
        leader = chosen[address ^ LEADER_BASE]
        return leader.bit_count() | ((leader & ((1 << PRIVATE_BITS) - 1)) << 8)
    if PARITY_BASE <= address < TAG_BASE:
        return LEADER_BASE | parity(address ^ PARITY_BASE)
    if TAG_BASE <= address < GAMMA_BASE:
        return displacement(address - TAG_BASE)
    if GAMMA_BASE <= address < PROFILE_BASE:
        return gamma(address - GAMMA_BASE)
    raise ValueError("not a fixed table address")


def operation_bound():
    extraction = 4 * (6 + 5)
    profile_build = 4 * SLOTS * 4
    pack_build = (1 + 3 * 8) * GROUPS * 25
    groups = GROUPS + GROUPS * SLOTS * 11
    improve = 125 * 41
    output = 7 + 4 * 6 + 6 * 4
    return BLOCKS * (extraction + profile_build + pack_build + groups + improve + output) + 16
