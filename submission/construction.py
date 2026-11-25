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
PARITY_BASE = 1 << 30
TAG_BASE = PARITY_BASE + (1 << PRIVATE_BITS)
SYNDROME_MASK = (1 << PARITY_BITS) - 1
TAG_MASK = SLOTS - 1
PRIVATE_MASK = (1 << PRIVATE_BITS) - 1
LOW_SPLIT = TAG_BASE + SLOTS
HIGH_SPLIT = LOW_SPLIT + (1 << 16)
LIST_BASE = HIGH_SPLIT + (1 << 15)
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
    return tuple(d | (weight << TAG_BITS) for weight, d in ranked[:RADIUS])


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


# The 12-bit field and the eight distinct evaluation coordinates.
MODULUS=0x1053
SIDES=8
DIMENSION=4
COMPONENT_BITS=SIDES*SIDE_BITS
PAYLOAD_BITS=SIDES*PRIVATE_BITS+DIMENSION*TAG_BITS
COMPONENTS=5
ANCHORS=tuple(combinations(range(SIDES),DIMENSION))
RADIUS=2
FOLLOW=(1<<PRIVATE_BITS)*RADIUS
TAG_SCRATCH=LIST_BASE+FOLLOW
TAG_SCRATCH_END=TAG_SCRATCH+SIDES*RADIUS
COEFF_BASE=1<<25
ROW_STRIDE=4*SLOTS
COEFF_END=COEFF_BASE+len(ANCHORS)*DIMENSION*ROW_STRIDE
DESCRIPTOR_BASE=1<<26
DESCRIPTOR_END=DESCRIPTOR_BASE+len(ANCHORS)*RADIUS**DIMENSION*DIMENSION
CACHED_RANKS=(2, 2, 2, 2)


def gf_mul(x,y):
    result=0
    while y:
        if y&1:
            result^=x
        y>>=1
        x<<=1
        if x&SLOTS:
            x^=MODULUS
    return result


def gf_pow(x,n):
    result=1
    while n:
        if n&1:
            result=gf_mul(result,x)
        x=gf_mul(x,x)
        n>>=1
    return result


EVALUATION_POINTS=tuple(x for x in range(SLOTS) if gf_pow(x,8)==x)


@lru_cache(None)
def coefficients(anchors):
    matrix=[]
    for i in anchors:
        denom=1
        for j in anchors:
            if j!=i:
                denom=gf_mul(denom,EVALUATION_POINTS[i]^EVALUATION_POINTS[j])
        inv=gf_pow(denom,SLOTS-2)
        assert gf_mul(denom,inv)==1
        row=[]
        for point in EVALUATION_POINTS:
            value=inv
            for j in anchors:
                if j!=i:
                    value=gf_mul(value,point^EVALUATION_POINTS[j])
            row.append(value)
        matrix.append(tuple(row))
    return tuple(matrix)


def tag_values(anchors, tags):
    result=[0]*SIDES
    for tag,row in zip(tags,coefficients(anchors),strict=True):
        for side,c in enumerate(row):
            result[side]^=gf_mul(c,tag)
    return tuple(result)


def encode(message):
    result=0
    for block in range(COMPONENTS):
        payload=(message>>(PAYLOAD_BITS*block))&((1<<PAYLOAD_BITS)-1)
        raw_tags=payload>>(SIDES*PRIVATE_BITS)
        tags=tag_values(ANCHORS[0],tuple((raw_tags>>(TAG_BITS*i))&TAG_MASK
                                       for i in range(DIMENSION)))
        for side,tag in enumerate(tags):
            u=(payload>>(PRIVATE_BITS*side))&PRIVATE_MASK
            word=u|((parity(u)^displacement(tag))<<PRIVATE_BITS)
            result|=word<<(COMPONENT_BITS*block+SIDE_BITS*side)
    return result


def coefficient_key(group,which):
    return COEFF_BASE+(group*DIMENSION+which)*ROW_STRIDE


def list_address(side,rank):
    return TAG_SCRATCH+side*RADIUS+rank


def descriptor_base(group, ranks):
    key=0
    for rank in ranks:
        key=key*RADIUS+rank
    return DESCRIPTOR_BASE+(group*RADIUS**DIMENSION+key)*DIMENSION


def table_word(address, chosen):
    if COEFF_BASE<=address<COEFF_END:
        row,offset=divmod(address-COEFF_BASE,ROW_STRIDE)
        group,i=divmod(row,DIMENSION)
        plane,tag=divmod(offset,SLOTS)
        anchors=ANCHORS[group]
        values=[gf_mul(c,tag) for c in coefficients(anchors)[i]]
        if plane<2:
            other=[s for s in range(SIDES) if s not in anchors]
            k,l=other[plane*2:plane*2+2]
            return displacement(values[k])|(displacement(values[l])<<32)
        return sum(values[s]<<(TAG_BITS*(s%4))
                   for s in range(4*(plane-2),4*(plane-1)))
    if DESCRIPTOR_BASE<=address<DESCRIPTOR_END:
        index,which=divmod(address-DESCRIPTOR_BASE,DIMENSION)
        group,digits=divmod(index,RADIUS**DIMENSION)
        rank=(digits//RADIUS**(DIMENSION-1-which))%RADIUS
        side=ANCHORS[group][which]
        plane2=coefficient_key(group,which)+2*SLOTS
        return list_address(side,rank)|(plane2<<32)
    if LEADER_BASE<=address<LEADER_BASE+(1<<PARITY_BITS):
        error=chosen[address^LEADER_BASE]
        return error.bit_count()|((error&PRIVATE_MASK)<<8)
    if PARITY_BASE<=address<TAG_BASE:
        return LEADER_BASE|parity(address^PARITY_BASE)
    if TAG_BASE<=address<LOW_SPLIT:
        return displacement(address-TAG_BASE)
    if LOW_SPLIT<=address<HIGH_SPLIT:
        return split_syndrome(address-LOW_SPLIT)
    if HIGH_SPLIT<=address<LIST_BASE:
        return split_syndrome((address-HIGH_SPLIT)<<16)
    if LIST_BASE<=address<TAG_SCRATCH:
        q,rank=divmod(address-LIST_BASE,RADIUS)
        return quotient_tags(q,chosen)[rank]
    raise ValueError("outside preloaded storage")


def side_list(block,side):
    yield from extract(COMPONENT_BITS*block+SIDE_BITS*side)
    yield "and_", f"U{side}", "Part", PRIVATE_MASK
    yield "shr", "High", "Part", PRIVATE_BITS
    yield "or_", "Addr", f"U{side}", PARITY_BASE
    yield "load", "Addr", "Addr"
    yield "xor", f"S{side}", "Addr", "High"
    yield "and_", "Low", f"S{side}", 65535
    yield "shr", "High", f"S{side}", 16
    yield "add", "Low", "Low", LOW_SPLIT
    yield "load", "Low", "Low"
    yield "add", "High", "High", HIGH_SPLIT-(1<<15)
    yield "load", "High", "High"
    yield "xor", "Split", "Low", "High"
    yield "and_", "Shift", "Split", TAG_MASK
    yield "shr", "Quot", "Split", TAG_BITS
    yield "mul", "List", "Quot", RADIUS
    yield "add", "List", "List", LIST_BASE
    for rank in range(RADIUS):
        if rank:
            yield "add", "Addr", "List", rank
            address="Addr"
        else:
            address="List"
        yield "load", "Entry", address
        yield "xor", "Entry", "Entry", "Shift"
        yield "store", None, list_address(side,rank), "Entry"


def contribution(group,which,side,rank,prefix):
    yield "load", "Entry", list_address(side,rank)
    yield "and_", "Key", "Entry", TAG_MASK
    yield "shr", prefix+"W", "Entry", TAG_BITS
    yield "or_", "Key", "Key", coefficient_key(group,which)
    yield "load", prefix+"L", "Key"
    yield "or_", "Key", "Key", SLOTS
    yield "load", prefix+"H", "Key"


def candidate_score(group,rank):
    for part in ("L","H"):
        yield "move", "Both", part+str(DIMENSION-1)
        yield "and_", "Addr", "Both", (1<<32)-1
        yield "load", "Entry", "Addr"
        yield "add", "Distance", "Distance", "Entry"
        yield "shr", "Addr", "Both", 32
        yield "load", "Entry", "Addr"
        yield "add", "Distance", "Distance", "Entry"
    yield "and_", "Distance", "Distance", 255
    key=descriptor_base(group,rank)
    if group==0 and not any(rank):
        yield "move", "Best", "Distance"
        yield "move", "Winner", key
    else:
        yield "ltu", "Better", "Distance", "Best"
        yield "select", "Best", "Better", "Distance", "Best"
        yield "select", "Winner", "Better", key, "Winner"


def candidate_tree(group,anchors,rank):
    which=len(rank)
    if which==DIMENSION:
        yield "move", "Distance", "W"+str(which-1)
        yield from candidate_score(group,rank)
        return
    for r in range(RADIUS):
        if r<CACHED_RANKS[which]:
            prefix=f"Cached{which}_{r}"
        else:
            prefix="Current"
            yield from contribution(group,which,anchors[which],r,prefix)
        for part in ("L","H","W"):
            dest=part+str(which)
            if which:
                yield ("add" if part=="W" else "xor"),dest,part+str(which-1),prefix+part
            elif part=="W":
                yield "move",dest,prefix+part
            else:
                yield "xor",dest,"Marked"+part,prefix+part
        yield from candidate_tree(group,anchors,rank+(r,))


def score_group(group,anchors):
    other=[side for side in range(SIDES) if side not in anchors]
    for pair,kind in enumerate(("L","H")):
        k,l=other[2*pair:2*pair+2]
        yield "shl","Marked"+kind,f"S{l}",32
        yield "or_","Marked"+kind,"Marked"+kind,f"S{k}"
    for which,count in enumerate(CACHED_RANKS):
        for rank in range(count):
            yield from contribution(group,which,anchors[which],rank,f"Cached{which}_{rank}")
    yield from candidate_tree(group,anchors,())


def winner_tags():
    yield "move","TagsL",0
    yield "move","TagsH",0
    for which in range(DIMENSION):
        yield "add","Addr","Winner",which
        yield "load","Descriptor","Addr"
        yield "and_","Addr","Descriptor",(1<<32)-1
        yield "shr","Key","Descriptor",32
        yield "load","Entry","Addr"
        yield "and_","Entry","Entry",TAG_MASK
        yield "or_","Key","Key","Entry"
        yield "load","Entry","Key"
        yield "xor","TagsL","TagsL","Entry"
        yield "or_","Key","Key",SLOTS
        yield "load","Entry","Key"
        yield "xor","TagsH","TagsH","Entry"


def component_schedule(block):
    start=block*PAYLOAD_BITS
    if start%64:
        yield "load",f"O{start//64}",32+start//64
    for side in range(SIDES):
        yield from side_list(block,side)
    for group,anchors in enumerate(ANCHORS):
        yield from score_group(group,anchors)
    yield from winner_tags()
    offset=block*PAYLOAD_BITS
    yield "move","OutputTags","TagsL"
    for side in range(SIDES):
        tags="TagsL" if side<4 else "TagsH"
        yield "and_","Tag",tags,TAG_MASK
        if side%4!=3:
            yield "shr",tags,tags,TAG_BITS
        yield "add","Addr","Tag",TAG_BASE
        yield "load","Addr","Addr"
        yield "xor","Addr","Addr",f"S{side}"
        yield "load","Entry","Addr"
        yield "shr","Error","Entry",8
        yield "xor","Decoded","Error",f"U{side}"
        yield from insert("Decoded",offset+PRIVATE_BITS*side,PRIVATE_BITS)
    for i in range(DIMENSION):
        yield "and_","Tag","OutputTags",TAG_MASK
        if i!=DIMENSION-1:
            yield "shr","OutputTags","OutputTags",TAG_BITS
        yield from insert("Tag",offset+SIDES*PRIVATE_BITS+TAG_BITS*i,TAG_BITS)
    stop=start+PAYLOAD_BITS
    if stop%64:
        yield "store",None,32+stop//64,f"O{stop//64}"
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
    for block in range(COMPONENTS):
        yield from allocate_registers(component_schedule(block))
