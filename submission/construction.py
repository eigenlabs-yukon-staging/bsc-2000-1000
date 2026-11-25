"""Fixed field, parity cells and decoder RAM tables."""
from qr import N, ORBITS, PERMUTATIONS, ROW_CAP, orbit_labels, parity as qr_parity, permute_byte, row_words

M = 125
ANCHOR_BITS = M
GROUPS = 8
WORDS = 2
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 53) | (1 << 54) | (1 << 57)
MODULUS = (1 << M) | MOD_LOW
RADIUS = 28
PRECISION = 16
BUDGETS = (2255, 2213, 2171, 2126, 2083, 2041, 1997, 1955, 1908, 1866, 1819, 1778, 1733, 1687, 1645, 1601, 1554, 1508, 1462, 1416, 1370, 1323, 1276, 1229, 1180, 1133, 1086, 1037, 987)
LIMIT = 2750
WIDTHS = (24,)*5 + (4,1)
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
MULT = FRAMES+6*TOTAL_CELLS
REDUCE_1 = MULT+16*WORDS
REDUCE_4 = REDUCE_1+2
SQUARE_BASE = REDUCE_4+16
PARITY_BASE = {4:SQUARE_BASE+BYTES*256*WORDS}
PARITY_BASE[1] = PARITY_BASE[4] + (1<<4)
LIST_BASE = {4:PARITY_BASE[1]+(1<<1)}
LIST_BASE[1] = LIST_BASE[4]+(1<<8)
QR_PARITY_BASE = LIST_BASE[1] + (1<<2)
QR_ORBIT_BASE = QR_PARITY_BASE + N
QR_PERM_BASE = QR_ORBIT_BASE + N
QR_ROW_BASE = QR_PERM_BASE + PERMUTATIONS * 6 * 256
BUDGET_BASE = QR_ROW_BASE + ORBITS * ROW_CAP
PROFILE_COUNTS = {24:14,4:3,1:2}
PROFILE_BASE = {24:BUDGET_BASE+len(BUDGETS)}
PROFILE_BASE[4] = PROFILE_BASE[24]+ORBITS
PROFILE_BASE[1] = PROFILE_BASE[4]+16
CAP_BASE = PROFILE_BASE[1]+2
CAP_OFFSETS = tuple(CAP_BASE+(LIMIT+1)*sum(PROFILE_COUNTS[w]
                    for w in (WIDTHS*(GROUPS-1))[:i]) for i in range(TOTAL_CELLS+1))
MEMORY_END = CAP_OFFSETS[-1]

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
        value=(value&MASK)^hi^(hi<<53)^(hi<<54)^(hi<<57)
    return value


def square(value):
    octets=value.to_bytes((value.bit_length()+7)//8,"little")
    return remainder(int.from_bytes(b"".join(SQUARE_BYTES[i] for i in octets),"little"))


def cube(value):
    return remainder(carryless(value,square(value)))


def hamming_parity(value):
    columns = (14,13,11,7)
    result = 0
    for i,column in enumerate(columns):
        if value >> i & 1:
            result ^= column
    return result


PARITY = {24:qr_parity,4:hamming_parity,1:int}


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
    if SQUARE_BASE<=address<PARITY_BASE[4]:
        row,word=divmod(address-SQUARE_BASE,WORDS)
        index,key=divmod(row,256)
        return (square_entry(index,key)>>(64*word))&MASK64
    for width in (4,1):
        offset = address-PARITY_BASE[width]
        if 0<=offset<1<<width:
            return PARITY[width](offset)
        offset = address-LIST_BASE[width]
        if 0<=offset<1<<(2*width):
            syndrome,rank=divmod(offset,1<<width)
            return candidate_row(width,syndrome)[rank]
    if QR_PARITY_BASE<=address<QR_ORBIT_BASE:
        return qr_parity(address-QR_PARITY_BASE)
    if QR_ORBIT_BASE<=address<QR_PERM_BASE:
        return orbit_labels()[0][address-QR_ORBIT_BASE]
    if QR_PERM_BASE<=address<QR_ROW_BASE:
        row,value=divmod(address-QR_PERM_BASE,256)
        perm,byte=divmod(row,6)
        return permute_byte(perm,byte,value)
    if QR_ROW_BASE<=address<BUDGET_BASE:
        orbit,ordinal=divmod(address-QR_ROW_BASE,ROW_CAP)
        return row_words(orbit,PRECISION)[ordinal]
    if BUDGET_BASE<=address<PROFILE_BASE[24]:
        return BUDGETS[address-BUDGET_BASE]
    if PROFILE_BASE[24]<=address<MEMORY_END:
        from dynamics import adaptive_word
        return adaptive_word(address)
    raise ValueError("not a fixed table address")
