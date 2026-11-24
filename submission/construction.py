"""Fixed linear code, table functions and primitive count."""
from collections import deque
from itertools import combinations

PRIVATE_BITS = 20
TAG_BITS = 10
PARITY_BITS = 30
SIDE_BITS = 50
BLOCKS = 20
WEIGHT_BASE = 1 << 30
CORRECTION_BASE = 2 << 30
PARITY_BASE = 3 << 30
TAG_BASE = PARITY_BASE + (1 << PRIVATE_BITS)
LIST_SIZE = 64
LIST_BASE = 1 << 26
LIST_END = LIST_BASE + (1 << PRIVATE_BITS) * LIST_SIZE

# Images of the 20 private basis bits, followed by the 10 tag basis bits.
COLUMNS = (222117643, 731568445, 1073074300, 951737589, 931311334, 14218408, 848738, 1046508680, 203304851, 378194437, 1046954136, 93887018, 727080981, 226724972, 373131137, 888004802, 324218201, 289108473, 928441295, 664932838, 172257281, 444849154, 453805060, 403846152, 506418192, 576033824, 78325824, 156651648, 1017617664, 121898496)


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


def encode(message):
    word = 0
    private_mask = (1 << PRIVATE_BITS) - 1
    for block in range(BLOCKS):
        payload = (message >> (50 * block)) & ((1 << 50) - 1)
        shift = displacement(payload >> 40)
        for side in range(2):
            u = (payload >> (PRIVATE_BITS * side)) & private_mask
            word |= (u | ((parity(u) ^ shift) << PRIVATE_BITS)) << (100 * block + SIDE_BITS * side)
    return word


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


def quotient(syndrome):
    tag_coordinate = syndrome & ((1 << TAG_BITS) - 1)
    return (syndrome ^ displacement(tag_coordinate)) >> TAG_BITS


def list_displacements(syndrome):
    columns = syndrome_columns()
    recorded = 0
    for weight in range(SIDE_BITS + 1):
        for coordinates in combinations(range(SIDE_BITS), weight):
            value = 0
            for i in coordinates:
                value ^= columns[i]
            if quotient(value) != syndrome:
                continue
            yield displacement(value & ((1 << TAG_BITS) - 1))
            recorded += 1
            if recorded == LIST_SIZE:
                return


def table_word(address, chosen):
    if LIST_BASE <= address < LIST_END:
        syndrome, slot = divmod(address - LIST_BASE, LIST_SIZE)
        for index, shift in enumerate(list_displacements(syndrome)):
            if index == slot:
                return shift
    if WEIGHT_BASE <= address < CORRECTION_BASE:
        return chosen[address ^ WEIGHT_BASE].bit_count()
    if CORRECTION_BASE <= address < PARITY_BASE:
        return chosen[address ^ CORRECTION_BASE] & ((1 << PRIVATE_BITS) - 1)
    if PARITY_BASE <= address < TAG_BASE:
        return WEIGHT_BASE | parity(address ^ PARITY_BASE)
    if TAG_BASE <= address < TAG_BASE + (1 << TAG_BITS):
        return displacement(address - TAG_BASE)
    raise ValueError("not a table address")


def operation_bound():
    extract = 2 * (6 + 5)
    list_search = 2 * (7 + LIST_SIZE * 9) + 61
    output = 1 + 2 * 4 + 4 + 4
    return BLOCKS * (extract + list_search + output) + 16
