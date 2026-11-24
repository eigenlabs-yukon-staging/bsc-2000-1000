"""Fixed field, parity cells and decoder RAM tables."""
M = 334
ANCHOR_BITS = 332
GROUPS = 3
WORDS = 6
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 1) | (1 << 2) | (1 << 5)
MODULUS = (1 << M) | MOD_LOW
RADIUS = 68
PRECISION = 16
EDGE_LIMIT = 1620
SLOPE = 40
LIMIT = EDGE_LIMIT + SLOPE*RADIUS
WIDTHS = (12,)*27 + (4,)*2 + (1,)*2
OFFSETS = tuple(sum(WIDTHS[:j]) for j in range(len(WIDTHS)))
SYS_WIDTHS = (ANCHOR_BITS,) + (M,)*(GROUPS-1)
PAR_WIDTHS = (M,ANCHOR_BITS) + (M,)*(GROUPS-2)
TOTAL_CELLS = (GROUPS - 1) * len(WIDTHS)
# A zero in this schedule denotes a width-one systematic-only slot.
SEARCH_WIDTHS = WIDTHS[:-2] + (0,0) + WIDTHS*(GROUPS-2)
SQUARE_BYTES = tuple(
    sum(((value >> bit) & 1) << (2*bit) for bit in range(8)).to_bytes(2,"little")
    for value in range(256)
)


ARRAYS = {name:48+WORDS*i for i,name in enumerate(
    tuple(f"{label}{j}" for label in ("Y","P","X","C","B") for j in range(GROUPS))
    + ("SQ","AP"))}
NEXT = 48 + WORDS*len(ARRAYS)
LAST = NEXT+M+1
SELECTED = LAST+M+1
FRAMES = SELECTED+TOTAL_CELLS
MULT = FRAMES+5*TOTAL_CELLS
REDUCE_1 = MULT+16*WORDS
REDUCE_4 = REDUCE_1+2
SQUARE_BASE = REDUCE_4+16
PARITY_BASE = {12:SQUARE_BASE+BYTES*256*WORDS}
PARITY_BASE[4] = PARITY_BASE[12] + (1<<12)
PARITY_BASE[1] = PARITY_BASE[4] + (1<<4)
LIST_BASE = {12:PARITY_BASE[1]+(1<<1)}
LIST_BASE[4] = LIST_BASE[12] + (1<<24)
LIST_BASE[1] = LIST_BASE[4] + (1<<8)
SYS_ONLY_BASE = LIST_BASE[1] + (1<<2)
MEMORY_END = SYS_ONLY_BASE + 2

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
        value=(value&MASK)^hi^(hi<<1)^(hi<<2)^(hi<<5)
    return value


def square(value):
    octets=value.to_bytes((value.bit_length()+7)//8,"little")
    return remainder(int.from_bytes(b"".join(SQUARE_BYTES[i] for i in octets),"little"))


def cube(value):
    return remainder(carryless(value,square(value)))


def golay_parity(value):
    generator = sum(1 << i for i in (0, 2, 4, 5, 6, 10, 11))
    remainder = poly_remainder(value << 11,generator)
    return remainder | (((value.bit_count() + remainder.bit_count()) % 2) << 11)


def hamming_parity(value):
    columns = (14,13,11,7)
    result = 0
    for i,column in enumerate(columns):
        if value >> i & 1:
            result ^= column
    return result


PARITY = {12:golay_parity,4:hamming_parity,1:int}


def linear(value):
    result = 0
    for w,i in zip(WIDTHS,OFFSETS,strict=True):
        result |= PARITY[w]((value>>i)&((1<<w)-1)) << i
    return result


def encode(message):
    sys_fields=[]
    offset=0
    for width in SYS_WIDTHS:
        sys_fields.append((message>>offset)&((1<<width)-1))
        offset+=width
    val=message
    for i,width in enumerate(PAR_WIDTHS):
        parity=linear(sys_fields[i])^cube(sys_fields[i-1])
        val|=(parity&((1<<width)-1))<<offset
        offset+=width
    return val


def cost(normalizer, weight):
    return ((normalizer**PRECISION-1)//weight**PRECISION).bit_length()


def square_entry(index,key):
    return square(key<<(8*index))


def candidate_row(width,syndrome):
    parity = PARITY[width]
    choices = sorted((e.bit_count()+(parity(e)^syndrome).bit_count(),e)
                     for e in range(1<<width))
    weights = [9**(2*width-d) for d in range(2*width+1)]
    normalizer = sum(weights[d] for d,_ in choices)
    costs = [cost(normalizer,weight) for weight in weights]
    return tuple(e|(d<<12)|(costs[d]<<17) for d,e in choices)


def table_word(address):
    if REDUCE_1<=address<REDUCE_4:
        return carryless(address-REDUCE_1,MOD_LOW)
    if REDUCE_4<=address<SQUARE_BASE:
        return carryless(address-REDUCE_4,MOD_LOW)
    if SQUARE_BASE<=address<PARITY_BASE[12]:
        row,word=divmod(address-SQUARE_BASE,WORDS)
        index,key=divmod(row,256)
        return (square_entry(index,key)>>(64*word))&MASK64
    for width in (12,4,1):
        offset = address-PARITY_BASE[width]
        if 0<=offset<1<<width:
            return PARITY[width](offset)
        offset = address-LIST_BASE[width]
        if 0<=offset<1<<(2*width):
            syndrome,rank=divmod(offset,1<<width)
            return candidate_row(width,syndrome)[rank]
    if SYS_ONLY_BASE<=address<MEMORY_END:
        e=address-SYS_ONLY_BASE
        return e|(e<<12)|(cost(10,9**(1-e))<<17)
    raise ValueError("not a fixed table address")
