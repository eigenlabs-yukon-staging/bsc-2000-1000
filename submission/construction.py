"""Four-field cyclic encoder, preloaded RAM tables, and exact count."""
from math import comb

M = 250
GROUPS = 4
WORDS = 4
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 32) | (1 << 40) | (1 << 59)
MODULUS = (1 << M) | MOD_LOW
BANDS = ((0, 16, 295), (17, 18, 287), (19, 20, 282), (21, 22, 275), (23, 24, 269), (25, 26, 263), (27, 28, 257), (29, 31, 249))
RADIUS = BANDS[-1][1]
SLICE = 18
SLICES = (M + SLICE - 1) // SLICE
ENTRIES = 1 << SLICE
AFFINE_STRIDE = SLICES * ENTRIES * WORDS
SQUARE_BYTES = tuple(
    sum(((value >> bit) & 1) << (2*bit) for bit in range(8)).to_bytes(2,"little")
    for value in range(256)
)

ARRAY_NAMES = tuple(
    f"{letter}{i}" for letter in ("Y","P","X","C","I","B") for i in range(GROUPS)
) + ("SQ",)
ARRAYS = {name:48+WORDS*i for i,name in enumerate(ARRAY_NAMES)}
TRAVERSAL_BASE = 48 + WORDS*len(ARRAY_NAMES)
TREE_COUNT = 1 + 2*(GROUPS-1)
NEXT = tuple(TRAVERSAL_BASE+2*(M+1)*i for i in range(TREE_COUNT))
LAST = tuple(n+M+1 for n in NEXT)
LIST_BASE = TRAVERSAL_BASE+2*TREE_COUNT*(M+1)
LISTS = tuple(LIST_BASE+M*i for i in range(2*(GROUPS-1)))
SAVED_BASE = LIST_BASE+2*(GROUPS-1)*M
SAVED = tuple(SAVED_BASE+(M+1)*WORDS*i for i in range(2*(GROUPS-1)))
MULT = SAVED_BASE+2*(GROUPS-1)*(M+1)*WORDS
REDUCE_1 = MULT+16*WORDS
REDUCE_4 = REDUCE_1+2
SQUARE_BASE = REDUCE_4+16
DEPTH_BASE = SQUARE_BASE+BYTES*256*WORDS
AFFINE_BASE = DEPTH_BASE+RADIUS+1
MEMORY_END = AFFINE_BASE+M*AFFINE_STRIDE


def carryless(a,b):
    result=0
    while b:
        bit=b&-b
        result ^= a<<(bit.bit_length()-1)
        b ^= bit
    return result


def poly_remainder(a,b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def remainder(value):
    while value>>M:
        hi=value>>M
        value=(value&MASK)^hi^(hi<<32)^(hi<<40)^(hi<<59)
    return value


def square(value):
    octets=value.to_bytes((value.bit_length()+7)//8,"little")
    return remainder(int.from_bytes(b"".join(SQUARE_BYTES[i] for i in octets),"little"))


def cube(value):
    return remainder(carryless(value,square(value)))


def encode(message):
    sys=tuple((message>>(M*i))&MASK for i in range(GROUPS))
    par=tuple(sys[i]^cube(sys[i-1]) for i in range(GROUPS))
    return sum(value<<(M*i) for i,value in enumerate(sys+par))


def square_entry(index,key):
    return square(key<<(8*index))


def affine_entry(direction,index,key):
    e=1<<direction
    result=cube(e) if index==0 else 0
    for bit in range(SLICE):
        j=index*SLICE+bit
        if j<M and (key>>bit)&1:
            x=1<<j
            result^=remainder((square(e)<<j)^(square(x)<<direction))
    return result


def table_word(address):
    if REDUCE_1<=address<REDUCE_4:
        return carryless(address-REDUCE_1,MOD_LOW)
    if REDUCE_4<=address<SQUARE_BASE:
        return carryless(address-REDUCE_4,MOD_LOW)
    if SQUARE_BASE<=address<DEPTH_BASE:
        row,word=divmod(address-SQUARE_BASE,WORDS)
        index,key=divmod(row,256)
        return (square_entry(index,key)>>(64*word))&MASK64
    if DEPTH_BASE<=address<AFFINE_BASE:
        depth=address-DEPTH_BASE
        return next(K for a,b,K in BANDS if a<=depth<=b)
    if AFFINE_BASE<=address<MEMORY_END:
        i,offset=divmod(address-AFFINE_BASE,AFFINE_STRIDE)
        entry,word=divmod(offset,WORDS)
        slice_index,key=divmod(entry,ENTRIES)
        return (affine_entry(i,slice_index,key)>>(64*word))&MASK64
    raise ValueError("not a fixed table address")


def operation_bound():
    k=WORDS
    shift=3*k+3
    square_ops=BYTES*(5+3*k)+k
    multiply_ops=3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops=square_ops+multiply_ops
    crosses=sum((i*SLICE)%64+SLICE>64 for i in range(SLICES))
    affine=2+2*k+SLICES*(6+3*k)+3*crosses
    save=2*(2+3*k)
    setup=14*M+9+4*k+1
    middle=44+46+3+affine+save
    leaf=44+46+11*GROUPS*k+2+affine+save
    per_orientation=2*k+1
    for lo,hi,K in BANDS:
        nodes=sum(comb(M,i) for i in range(lo,hi+1))
        per_anchor=cube_ops+3+44
        visits=1
        for level in range(1,GROUPS):
            node_cost=leaf if level==GROUPS-1 else middle
            bound=1<<(K-6*(GROUPS-1-level))
            per_anchor += setup*visits + node_cost*bound
            visits=bound
        per_orientation += nodes*per_anchor
    total=100000+GROUPS*(cube_ops+per_orientation)
    return total,(square_ops,multiply_ops,affine,save,setup,middle,leaf,per_orientation)
