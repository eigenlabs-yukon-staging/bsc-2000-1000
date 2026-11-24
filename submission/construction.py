"""Fixed linear code, table functions and primitive count."""
from collections import deque

PRIVATE_BITS = 20
TAG_BITS = 10
PARITY_BITS = 30
SIDE_BITS = 50
BLOCKS = 10
GROUPS = 1 << (TAG_BITS - 3)
SLOTS = 1 << TAG_BITS
WEIGHT_BASE = 1 << 30
CORRECTION_BASE = 2 << 30
PARITY_BASE = 3 << 30
TAG_BASE = PARITY_BASE + (1 << PRIVATE_BITS)
GAMMA_BASE = TAG_BASE + SLOTS
PROFILE_BASE = GAMMA_BASE + SLOTS
PACKED_BASE = PROFILE_BASE + 4 * SLOTS
END = PACKED_BASE + (1 + 3 * 8) * GROUPS
HIGH = 0x8080808080808080
REPLICATE = 0x0101010101010101

# Images of the 20 private basis bits, followed by the 10 tag basis bits.
COLUMNS = (591308903, 413046786, 804950232, 668036205, 794244845, 293736808, 97193133, 193452214, 315598187, 437487036, 410678129, 1002135223, 190165052, 775971152, 405512570, 119709043, 817603403, 89359809, 1055644720, 100447513, 709128193, 713284610, 453805060, 672281608, 1043289104, 844469280, 78325824, 693522560, 212311296, 658769408)


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
    result = 0
    multiplier = 237
    while multiplier:
        if multiplier & 1:
            result ^= tag
        tag = ((tag << 1) & (SLOTS - 1)) ^ (9 if tag >> 9 else 0)
        multiplier >>= 1
    return result


def tag_tuple(a, s):
    d = gamma(s)
    return a, a ^ s, a ^ d, a ^ s ^ d


def encode(message):
    word = 0
    private_mask = (1 << PRIVATE_BITS) - 1
    for block in range(BLOCKS):
        payload = (message >> (100 * block)) & ((1 << 100) - 1)
        tags = tag_tuple((payload >> 80) & (SLOTS - 1), payload >> 90)
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
    if WEIGHT_BASE <= address < CORRECTION_BASE:
        return chosen[address ^ WEIGHT_BASE].bit_count()
    if CORRECTION_BASE <= address < PARITY_BASE:
        return chosen[address ^ CORRECTION_BASE] & ((1 << PRIVATE_BITS) - 1)
    if PARITY_BASE <= address < TAG_BASE:
        return WEIGHT_BASE | parity(address ^ PARITY_BASE)
    if TAG_BASE <= address < TAG_BASE + (1 << TAG_BITS):
        return WEIGHT_BASE ^ CORRECTION_BASE ^ displacement(address - TAG_BASE)
    if GAMMA_BASE <= address < PROFILE_BASE:
        return gamma(address - GAMMA_BASE)
    raise ValueError("not a fixed table address")


def operation_bound():
    extraction = 4 * (6 + 5)
    profile_build = 4 * SLOTS * 3
    pack_build = (1 + 3 * 8) * GROUPS * 25
    groups = GROUPS + GROUPS * SLOTS * 11
    improve = 121 * 41
    output = 7 + 4 * 5 + 6 * 4
    return BLOCKS * (extraction + profile_build + pack_build + groups + improve + output) + 16
