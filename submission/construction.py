"""Fixed linear code, table functions and primitive count."""
from collections import deque
from functools import lru_cache
from itertools import combinations

PRIVATE_BITS = 19
TAG_BITS = 12
PARITY_BITS = 31
SIDE_BITS = 50
BLOCKS = 10
SLOTS = 1 << TAG_BITS
LEADER_BASE = 1 << 31
PAIR_BASE = 1 << 27
PAIR_PLANE = 1 << 24
PAIR_STRIDE = 4 * PAIR_PLANE
PAIRS = tuple(combinations(range(4),2))
PAIR_END = PAIR_BASE + len(PAIRS) * PAIR_STRIDE
PARITY_BASE = 1 << 30
TAG_BASE = PARITY_BASE + (1 << PRIVATE_BITS)
GAMMA_BASE = TAG_BASE + SLOTS
SIDE_LIST_COUNTS = (1, 2, 2, 2)
LIST_COUNT = max(SIDE_LIST_COUNTS)
SYNDROME_MASK = (1 << PARITY_BITS) - 1
TAG_MASK = SLOTS - 1
PRIVATE_MASK = (1 << PRIVATE_BITS) - 1
RECORD_BITS = 18
LIST_SHIFT = 8
CORRECTION_SHIFT = 45


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


GAMMA_COLUMNS = (2214, 3535, 2322, 2924, 1333, 2150, 1796, 2021, 533, 1940, 3031, 1017)


def gamma(tag):
    return linear_map(tag, GAMMA_COLUMNS)


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


@lru_cache(None)
def split_basis():
    selected = []
    pivots = {}
    for original in (*COLUMNS[PRIVATE_BITS:], *(1 << i for i in range(PARITY_BITS))):
        reduced = original
        mark = 1 << len(selected)
        while reduced:
            bit = reduced.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = reduced, mark
                selected.append(original)
                break
            vector, prior = pivots[bit]
            reduced ^= vector
            mark ^= prior
    assert len(selected) == PARITY_BITS
    assert tuple(selected[:TAG_BITS]) == COLUMNS[PRIVATE_BITS:]
    return tuple(selected), pivots


def split_syndrome(syndrome):
    _, pivots = split_basis()
    result = 0
    while syndrome:
        vector, mark = pivots[syndrome.bit_length() - 1]
        syndrome ^= vector
        result ^= mark
    return result


def quotient_tags(quotient, chosen):
    basis, _ = split_basis()
    representative = linear_map(quotient, basis[TAG_BITS:])
    ranked = sorted((chosen[representative ^ displacement(d)].bit_count(), d)
                    for d in range(SLOTS))
    return tuple(d | (weight << TAG_BITS) for weight, d in ranked[:LIST_COUNT])


def table_word(address, chosen):
    if PAIR_BASE <= address < PAIR_END:
        pair_id, offset = divmod(address-PAIR_BASE,PAIR_STRIDE)
        plane, key = divmod(offset,PAIR_PLANE)
        i,j = PAIRS[pair_id]
        u,v = divmod(key,SLOTS)
        k,l = pair_slots(i,j)
        d = gamma(u^v)
        tags = [0]*4
        tags[i],tags[j],tags[k],tags[l] = u,v,u^d,v^d
        if plane == 0:
            return displacement(tags[k]) | (displacement(tags[l]) << 32)
        if plane == 2:
            return sum(tag << (TAG_BITS*side) for side,tag in enumerate(tags))
        return 0
    if LEADER_BASE <= address < LEADER_BASE + (1 << PARITY_BITS):
        syndrome = address ^ LEADER_BASE
        leader = chosen[syndrome]
        coordinates = split_syndrome(syndrome)
        q,shift = divmod(coordinates,SLOTS)
        word = leader.bit_count() | ((leader & PRIVATE_MASK) << CORRECTION_SHIFT)
        for rank,relative in enumerate(quotient_tags(q,chosen)):
            absolute = relative ^ shift
            word |= absolute << (LIST_SHIFT+RECORD_BITS*rank)
        return word
    if PARITY_BASE <= address < TAG_BASE:
        return LEADER_BASE | parity(address ^ PARITY_BASE)
    if TAG_BASE <= address < GAMMA_BASE:
        return displacement(address - TAG_BASE)
    raise ValueError("not a preloaded address")


def pair_slots(i, j):
    # For these two known tags, slot k has ti xor G(ti xor tj).
    test_tuple = tag_tuple(0, 1)
    first = test_tuple[i] ^ gamma(test_tuple[i] ^ test_tuple[j])
    second = test_tuple[j] ^ gamma(test_tuple[i] ^ test_tuple[j])
    return test_tuple.index(first), test_tuple.index(second)


def extract(start):
    word, shift = divmod(start, 64)
    yield "load", "Part", word
    if shift:
        yield "shr", "Part", "Part", shift
    if shift + SIDE_BITS > 64:
        yield "load", "Temp", word + 1
        yield "shl", "Temp", "Temp", 64 - shift
        yield "or_", "Part", "Part", "Temp"
    yield "and_", "Part", "Part", (1 << SIDE_BITS) - 1


def insert(value, start, width):
    word, shift = divmod(start, 64)
    out = f"O{word}"
    if shift:
        yield "shl", "Temp", value, shift
        yield "or_", out, out, "Temp"
    else:
        yield "move", out, value
    if shift + width >= 64:
        yield "store", None, 32 + word, out
    if shift + width > 64:
        next_out = f"O{word+1}"
        yield "shr", next_out, value, 64 - shift


def relative_list(block, side):
    yield from extract(200 * block + SIDE_BITS * side)
    yield "and_", f"U{side}", "Part", PRIVATE_MASK
    yield "shr", "High", "Part", PRIVATE_BITS
    yield "or_", "Addr", f"U{side}", PARITY_BASE
    yield "load", "Addr", "Addr"
    yield "xor", f"S{side}", "Addr", "High"
    yield "load", "List", f"S{side}"
    for rank in range(SIDE_LIST_COUNTS[side]):
        yield "shr", "Record", "List", LIST_SHIFT + RECORD_BITS * rank
        yield "and_", tag_register(side,rank), "Record", TAG_MASK
        yield "shr", "Weight", "Record", TAG_BITS
        yield "and_", distance_register(side,rank), "Weight", 63


def tag_register(side, rank):
    return f"TagList{side}_{rank}"


def distance_register(side, rank):
    return f"DistList{side}_{rank}"


def leader_entry(side, tag):
    yield "add", "Addr", tag, TAG_BASE
    yield "load", "Addr", "Addr"
    yield "xor", "Addr", f"S{side}", "Addr"
    yield "load", "Entry", "Addr"


def try_pairs(pair_id, i, j):
    k,l = pair_slots(i,j)
    base=PAIR_BASE+pair_id*PAIR_STRIDE
    yield "shl", "Syndromes", f"S{l}", 32
    yield "or_", "Syndromes", "Syndromes", f"S{k}"
    for u in range(SIDE_LIST_COUNTS[i]):
        first_tag=tag_register(i,u)
        first_dist=distance_register(i,u)
        yield "shl", "Prefix", first_tag, TAG_BITS
        yield "or_", "Prefix", "Prefix", base
        for v in range(SIDE_LIST_COUNTS[j]):
            second_tag=tag_register(j,v)
            second_dist=distance_register(j,v)
            yield "or_", "Key", "Prefix", second_tag
            yield "add", "Distance", first_dist, second_dist
            yield "load", "Both", "Key"
            yield "xor", "Both", "Both", "Syndromes"
            yield "and_", "Addr", "Both", (1 << 32) - 1
            yield "load", "Entry", "Addr"
            yield "add", "Distance", "Distance", "Entry"
            yield "shr", "Addr", "Both", 32
            yield "load", "Entry", "Addr"
            yield "add", "Distance", "Distance", "Entry"
            yield "and_", "Distance", "Distance", 255
            if (pair_id, u, v) == (0, 0, 0):
                yield "move", "Best", "Distance"
                yield "move", "Winner", "Key"
            else:
                yield "ltu", "Improve", "Distance", "Best"
                yield "select", "Best", "Improve", "Distance", "Best"
                yield "select", "Winner", "Improve", "Key", "Winner"


def virtual_schedule():
    for block in range(BLOCKS):
        for side in range(4):
            yield from relative_list(block, side)
        for pair_id,(i,j) in enumerate(PAIRS):
            yield from try_pairs(pair_id,i,j)
        yield "or_", "Winner", "Winner", 2*PAIR_PLANE
        yield "load", "Tags", "Winner"
        for side in range(4):
            yield "and_", f"B{side}", "Tags", TAG_MASK
            if side<3:
                yield "shr", "Tags", "Tags", TAG_BITS
        yield "xor", "Direction", "B0", "B1"
        for side in range(4):
            yield from leader_entry(side, f"B{side}")
            yield "shr", "Error", "Entry", CORRECTION_SHIFT
            yield "xor", "Decoded", f"U{side}", "Error"
            yield from insert("Decoded", 100 * block + PRIVATE_BITS * side, PRIVATE_BITS)
        yield from insert("B0", 100 * block + 76, TAG_BITS)
        yield from insert("Direction", 100 * block + 88, TAG_BITS)
    yield "store", None, 47, "O15"


def allocate_registers(instructions):
    from heapq import heapify, heappop, heappush

    history=[]
    latest={}
    last_use=[]
    for index,(op,dest,*args) in enumerate(instructions):
        sources=tuple(latest[arg] if isinstance(arg,str) else None for arg in args)
        for source in sources:
            if source is not None:
                last_use[source]=index
        if dest is None:
            result=None
        else:
            result=len(last_use)
            latest[dest]=result
            last_use.append(index)
        history.append((op,result,args,sources))

    free=list(range(64))
    heapify(free)
    live={}
    for index,(op,result,args,sources) in enumerate(history):
        operands=tuple(arg if source is None else f"R{live[source]}"
                       for arg,source in zip(args,sources,strict=True))
        for source in set(sources)-{None}:
            if last_use[source]==index:
                heappush(free,live.pop(source))
        dest=None
        if result is not None:
            slot=heappop(free)
            live[result]=slot
            dest=f"R{slot}"
        yield op,dest,*operands
        if result is not None and last_use[result]==index:
            heappush(free,live.pop(result))


def schedule():
    yield from allocate_registers(virtual_schedule())
