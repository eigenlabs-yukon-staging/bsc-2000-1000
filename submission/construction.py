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
COLUMNS = (1051358863, 590810211, 1475472909, 1978738269, 651388888, 232409661, 1369593774, 2076077886, 429576749, 416644663, 1299312799, 1275935853, 711053259, 784118576, 130579930, 548027298, 790821889, 618131721, 451792847, 344514561, 689029122, 705654788, 1815220232, 541642768, 2025672736, 1230393408, 1387044992, 1700348160, 849244672, 1561334784, 2048927744)


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


MODULUS=0x1053
SIDES=10
DIMENSION=5
COMPONENTS=4
COMPONENT_BITS=SIDES*SIDE_BITS
PAYLOAD_BITS=SIDES*PRIVATE_BITS+DIMENSION*TAG_BITS
EVALUATION_POINTS=tuple(range(SIDES))
ANCHORS=tuple(combinations(range(SIDES),DIMENSION))
RANK_COSTS=(0, 1, 2, 3, 4)
RADIUS=len(RANK_COSTS)
BUDGET=4
PAIRED=(DIMENSION+1)//2
CACHED_DEPTHS=1
UNCACHED_DEPTHS=DIMENSION-CACHED_DEPTHS

FIELD_DISP_BASE=1<<24
FIELD_TAG_BASE=1<<25
GROUP_BASE=1<<26
GROUP_STRIDE=DIMENSION*(SIDES+1)
GROUP_END=GROUP_BASE+len(ANCHORS)*GROUP_STRIDE

LIST_END=LIST_BASE+(1<<PRIVATE_BITS)*RADIUS
U_WORK=LIST_END
S_WORK=U_WORK+SIDES
SHORTLIST=S_WORK+SIDES
WEIGHT_WORK=SHORTLIST+SIDES*RADIUS
DISP_WORK=WEIGHT_WORK+UNCACHED_DEPTHS*RADIUS
BEST_SCORE=DISP_WORK+UNCACHED_DEPTHS*RADIUS*PAIRED
BEST_KEY=BEST_SCORE+1
WORK_END=BEST_KEY+1


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


@lru_cache(64)
def coefficients(anchors):
    matrix=[]
    for i in anchors:
        denominator=1
        for j in anchors:
            if i!=j:
                denominator=gf_mul(denominator,EVALUATION_POINTS[i]^EVALUATION_POINTS[j])
        inverse=gf_pow(denominator,SLOTS-2)
        assert gf_mul(denominator,inverse)==1
        row=[]
        for point in EVALUATION_POINTS:
            c=inverse
            for j in anchors:
                if j!=i:
                    c=gf_mul(c,point^EVALUATION_POINTS[j])
            row.append(c)
        matrix.append(tuple(row))
    return tuple(matrix)


def tag_values(anchors,tags):
    result=[0]*SIDES
    for tag,row in zip(tags,coefficients(anchors),strict=True):
        for side,c in enumerate(row):
            result[side]^=gf_mul(c,tag)
    return tuple(result)


def encode(message):
    result=0
    for block in range(COMPONENTS):
        payload=(message>>(PAYLOAD_BITS*block))&((1<<PAYLOAD_BITS)-1)
        raw=payload>>(SIDES*PRIVATE_BITS)
        tags=tag_values(ANCHORS[0],tuple((raw>>(TAG_BITS*i))&TAG_MASK
                                       for i in range(DIMENSION)))
        for side,tag in enumerate(tags):
            u=(payload>>(PRIVATE_BITS*side))&PRIVATE_MASK
            word=u|((parity(u)^displacement(tag))<<PRIVATE_BITS)
            result|=word<<(COMPONENT_BITS*block+SIDE_BITS*side)
    return result


def table_word(address,chosen):
    if FIELD_DISP_BASE<=address<FIELD_TAG_BASE+(1<<24):
        bank,index=divmod(address-FIELD_DISP_BASE,1<<24)
        c,t=divmod(index,SLOTS)
        result=gf_mul(c,t)
        return displacement(result) if bank==0 else result
    if GROUP_BASE<=address<GROUP_END:
        group,within=divmod(address-GROUP_BASE,GROUP_STRIDE)
        which,slot=divmod(within,SIDES+1)
        anchors=ANCHORS[group]
        if slot==0:
            return SHORTLIST+RADIUS*anchors[which]
        return FIELD_TAG_BASE+SLOTS*coefficients(anchors)[which][slot-1]
    if LEADER_BASE<=address<LEADER_BASE+(1<<PARITY_BITS):
        leader=chosen[address^LEADER_BASE]
        return leader.bit_count()|((leader&PRIVATE_MASK)<<8)
    if PARITY_BASE<=address<TAG_BASE:
        return LEADER_BASE|parity(address^PARITY_BASE)
    if TAG_BASE<=address<LOW_SPLIT:
        return displacement(address-TAG_BASE)
    if LOW_SPLIT<=address<HIGH_SPLIT:
        return split_syndrome(address-LOW_SPLIT)
    if HIGH_SPLIT<=address<LIST_BASE:
        return split_syndrome((address-HIGH_SPLIT)<<16)
    if LIST_BASE<=address<LIST_END:
        q,rank=divmod(address-LIST_BASE,RADIUS)
        return quotient_tags(q,chosen)[rank]
    raise ValueError("not a preloaded address")


def side_list(block,side):
    yield from extract(COMPONENT_BITS*block+SIDE_BITS*side)
    yield "and_", "U", "Part", PRIVATE_MASK
    yield "store", None, U_WORK+side, "U"
    yield "shr", "High", "Part", PRIVATE_BITS
    yield "or_", "Addr", "U", PARITY_BASE
    yield "load", "Addr", "Addr"
    yield "xor", "S", "Addr", "High"
    yield "store", None, S_WORK+side, "S"
    yield "and_", "Low", "S", 65535
    yield "shr", "High", "S", 16
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
        yield "store", None, SHORTLIST+side*RADIUS+rank, "Entry"


def setup_component(block):
    for side in range(SIDES):
        yield from side_list(block,side)


def weight_index(which,rank):
    return WEIGHT_WORK+which*RADIUS+rank


def displacement_index(which,rank,pair):
    return DISP_WORK+(which*RADIUS+rank)*PAIRED+pair


def prepare_group(anchors):
    missing=tuple(s for s in range(SIDES) if s not in anchors)
    for pair in range(PAIRED):
        k=missing[2*pair]
        yield "load",f"Marked{pair}",S_WORK+k
        if 2*pair+1<DIMENSION:
            l=missing[2*pair+1]
            yield "load","Hi",S_WORK+l
            yield "shl","Hi","Hi",32
            yield "or_",f"Marked{pair}",f"Marked{pair}","Hi"
    for which,(side,row) in enumerate(zip(anchors,coefficients(anchors),strict=True)):
        for rank in range(RADIUS):
            yield "load","Entry",SHORTLIST+side*RADIUS+rank
            yield "and_","Tag","Entry",TAG_MASK
            if which<UNCACHED_DEPTHS:
                yield "shr","Weight","Entry",TAG_BITS
                yield "store",None,weight_index(which,rank),"Weight"
            else:
                yield "shr",cached_weight(which,rank),"Entry",TAG_BITS
            for pair in range(PAIRED):
                k=missing[2*pair]
                yield "or_","Addr","Tag",FIELD_DISP_BASE+SLOTS*row[k]
                yield "load","Both","Addr"
                if 2*pair+1<DIMENSION:
                    l=missing[2*pair+1]
                    yield "or_","Addr","Tag",FIELD_DISP_BASE+SLOTS*row[l]
                    yield "load","Hi","Addr"
                    yield "shl","Hi","Hi",32
                    yield "or_","Both","Both","Hi"
                if which<UNCACHED_DEPTHS:
                    yield "store",None,displacement_index(which,rank,pair),"Both"
                else:
                    yield "move",cached_displacement(which,rank,pair),"Both"


def cached_weight(depth,rank):
    return f"CachedWeight{depth}_{rank}"


def cached_displacement(depth,rank,pair):
    return f"CachedDisp{depth}_{rank}_{pair}"


def score_leaf(group,ranks):
    yield "move","Distance",f"PrefixW{DIMENSION-1}"
    for pair in range(PAIRED):
        both=f"Prefix{DIMENSION-1}_{pair}"
        if 2*pair+1<DIMENSION:
            yield "and_","Addr",both,(1<<32)-1
            yield "load","Entry","Addr"
            yield "add","Distance","Distance","Entry"
            yield "shr","Addr",both,32
            address="Addr"
        else:
            address=both
        yield "load","Entry",address
        yield "add","Distance","Distance","Entry"
    yield "and_","Distance","Distance",255
    ordinal=0
    for rank in ranks:
        ordinal=ordinal*RADIUS+rank
    key=group*RADIUS**DIMENSION+ordinal
    if key==0:
        yield "move","Best","Distance"
        yield "move","Winner",0
    else:
        yield "ltu","Better","Distance","Best"
        yield "select","Best","Better","Distance","Best"
        yield "select","Winner","Better",key,"Winner"


def score_tree(group,ranks,budget):
    depth=len(ranks)
    if depth==DIMENSION:
        yield from score_leaf(group,ranks)
        return
    for rank,cost in enumerate(RANK_COSTS):
        if cost>budget:
            break
        if depth<UNCACHED_DEPTHS:
            local_weight="Weight"
            yield "load",local_weight,weight_index(depth,rank)
        else:
            local_weight=cached_weight(depth,rank)
        weight=f"PrefixW{depth}"
        if depth:
            yield "add",weight,f"PrefixW{depth-1}",local_weight
        else:
            yield "move",weight,local_weight
        for pair in range(PAIRED):
            if depth<UNCACHED_DEPTHS:
                displacement="Both"
                yield "load",displacement,displacement_index(depth,rank,pair)
            else:
                displacement=cached_displacement(depth,rank,pair)
            prefix=f"Prefix{depth}_{pair}"
            previous=f"Prefix{depth-1}_{pair}" if depth else f"Marked{pair}"
            yield "xor",prefix,previous,displacement
        yield from score_tree(group,ranks+(rank,),budget-cost)


def group_schedule(group,anchors):
    yield from prepare_group(anchors)
    if group:
        yield "load","Best",BEST_SCORE
        yield "load","Winner",BEST_KEY
    yield from score_tree(group,(),BUDGET)
    yield "store",None,BEST_SCORE,"Best"
    yield "store",None,BEST_KEY,"Winner"


def decode_winner(block):
    start=block*PAYLOAD_BITS
    if start%64:
        yield "load",f"O{start//64}",32+start//64
    yield "load","Winner",BEST_KEY
    yield "divu","Group","Winner",RADIUS**DIMENSION
    yield "mul","Group","Group",GROUP_STRIDE
    yield "add","Group","Group",GROUP_BASE
    for which in range(DIMENSION):
        yield "add","Addr","Group",which*(SIDES+1)
        yield "load","List","Addr"
        yield "divu","Rank","Winner",RADIUS**(DIMENSION-1-which)
        yield "modu","Rank","Rank",RADIUS
        yield "add","Addr","List","Rank"
        yield "load","Entry","Addr"
        yield "and_",f"Selected{which}","Entry",TAG_MASK
    for side in range(SIDES):
        tag=f"OutputTag{side}" if side<DIMENSION else "ScratchTag"
        for which in range(DIMENSION):
            yield "add","Addr","Group",which*(SIDES+1)+side+1
            yield "load","Key","Addr"
            yield "or_","Key","Key",f"Selected{which}"
            yield "load","Partial","Key"
            if which:
                yield "xor",tag,tag,"Partial"
            else:
                yield "move",tag,"Partial"
        yield "load","S",S_WORK+side
        yield "add","Addr",tag,TAG_BASE
        yield "load","Addr","Addr"
        yield "xor","Addr","S","Addr"
        yield "load","Entry","Addr"
        yield "shr","Error","Entry",8
        yield "load","U",U_WORK+side
        yield "xor","Decoded","U","Error"
        yield from insert("Decoded",start+side*PRIVATE_BITS,PRIVATE_BITS)
    for side in range(DIMENSION):
        yield from insert(f"OutputTag{side}",start+SIDES*PRIVATE_BITS+TAG_BITS*side,TAG_BITS)
    end=start+PAYLOAD_BITS
    if end%64:
        yield "store",None,32+end//64,f"O{end//64}"


def schedule():
    for block in range(COMPONENTS):
        yield from allocate_registers(setup_component(block))
        for group,anchors in enumerate(ANCHORS):
            yield from allocate_registers(group_schedule(group,anchors))
        yield from allocate_registers(decode_winner(block))


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


