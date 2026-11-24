"""Fixed field, parity cells and decoder RAM tables."""
from math import comb

M = 200
GROUPS = 5
WORDS = 4
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 11) | (1 << 19) | (1 << 30)
MODULUS = (1 << M) | MOD_LOW
RADIUS = 33
PRECISION = 8
LIMIT = 1144
WIDTHS = (12,)*16 + (4,)*2
OFFSETS = tuple(sum(WIDTHS[:j]) for j in range(len(WIDTHS)))
TOTAL_CELLS = (GROUPS - 1) * len(WIDTHS)
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
LIST_BASE = {12:PARITY_BASE[4]+(1<<4)}
LIST_BASE[4] = LIST_BASE[12] + (1<<24)
MEMORY_END = LIST_BASE[4] + (1<<8)

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
        value=(value&MASK)^hi^(hi<<11)^(hi<<19)^(hi<<30)
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


PARITY = {12:golay_parity,4:hamming_parity}


def linear(value):
    result = 0
    for w,i in zip(WIDTHS,OFFSETS,strict=True):
        result |= PARITY[w]((value>>i)&((1<<w)-1)) << i
    return result


def encode(message):
    sys=tuple((message>>(M*i))&MASK for i in range(GROUPS))
    par=tuple(linear(sys[i])^cube(sys[i-1]) for i in range(GROUPS))
    return sum(value<<(M*i) for i,value in enumerate(sys+par))


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
    for width in (12,4):
        offset = address-PARITY_BASE[width]
        if 0<=offset<1<<width:
            return PARITY[width](offset)
        offset = address-LIST_BASE[width]
        if 0<=offset<1<<(2*width):
            syndrome,rank=divmod(offset,1<<width)
            return candidate_row(width,syndrome)[rank]
    raise ValueError("not a fixed table address")


def operation_bound():
    k = WORDS
    shift = 3*k+3
    square_ops = BYTES*(5+3*k)+k
    multiply_ops = 3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops = square_ops+multiply_ops
    # Pseudocode gives uniform rounded-up costs per event.
    # Intermediate cell depths become complete fields at multiples of 18.
    counts = [1]+[1<<((LIMIT-(TOTAL_CELLS-ell)+PRECISION-1)//PRECISION)
                  for ell in range(1,TOTAL_CELLS+1)]
    per_anchor = 44+250+cube_ops
    for ell in range(1,TOTAL_CELLS+1):
        per_anchor += 40*counts[ell-1]+28*counts[ell]
        if ell%len(WIDTHS)==0:
            additional = 180 if ell<TOTAL_CELLS else 200
            per_anchor += (cube_ops+additional)*counts[ell]
    total = 1000+sum(comb(M,i) for i in range(RADIUS+1))*per_anchor
    return total
