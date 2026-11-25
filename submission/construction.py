"""Fixed field, parity cells and decoder RAM tables."""
from qr import N, ORBITS, PERMUTATIONS, ROW_CAP, orbit_labels, parity as qr_parity, permute_byte, row_words

M = 500
ANCHOR_BITS = M
GROUPS = 2
WORDS = 8
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 8) | (1 << 13) | (1 << 21)
MODULUS = (1 << M) | MOD_LOW
RADIUS = 102
PRECISION = 16
BUDGETS = (3865, 3865, 3865, 3865, 3865, 3865, 3865, 3865, 3865, 3734, 3733, 3733, 3728, 3710, 3710, 3710, 3690, 3671, 3671, 3671, 3671, 3671, 3671, 3628, 3628, 3570, 3526, 3505, 3467, 3430, 3391, 3362, 3300, 3264, 3226, 3167, 3149, 3098, 3061, 3018, 2965, 2924, 2894, 2857, 2811, 2759, 2727, 2691, 2653, 2611, 2563, 2527, 2490, 2451, 2411, 2370, 2330, 2293, 2250, 2215, 2177, 2139, 2101, 2064, 2021, 1981, 1942, 1905, 1867, 1824, 1788, 1745, 1706, 1666, 1627, 1587, 1546, 1504, 1465, 1424, 1383, 1341, 1300, 1261, 1219, 1177, 1136, 1094, 1053, 1010, 968, 925, 882, 840, 797, 754, 709, 666, 620, 578, 532, 487, 442)
LIMIT = 5600
WIDTHS = (24,)*20 + (12,4,4)
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
PARITY_BASE = {12:SQUARE_BASE+BYTES*256*WORDS}
PARITY_BASE[4] = PARITY_BASE[12] + (1<<12)
LIST_BASE = {12:PARITY_BASE[4]+(1<<4)}
LIST_BASE[4] = LIST_BASE[12] + (1<<24)
QR_PARITY_BASE = LIST_BASE[4] + (1<<8)
QR_ORBIT_BASE = QR_PARITY_BASE + N
QR_PERM_BASE = QR_ORBIT_BASE + N
QR_ROW_BASE = QR_PERM_BASE + PERMUTATIONS * 6 * 256
BUDGET_BASE = QR_ROW_BASE + ORBITS * ROW_CAP
PROFILE_COUNTS = {24:14,12:5,4:3}
PROFILE_BASE = {24:BUDGET_BASE+len(BUDGETS)}
PROFILE_BASE[12] = PROFILE_BASE[24]+ORBITS
PROFILE_BASE[4] = PROFILE_BASE[12]+4096
CAP_BASE = PROFILE_BASE[4]+16
CAP_OFFSETS = tuple(CAP_BASE+(LIMIT+1)*sum(PROFILE_COUNTS[w]
                    for w in WIDTHS[:i]) for i in range(TOTAL_CELLS+1))
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
        value=(value&MASK)^hi^(hi<<8)^(hi<<13)^(hi<<21)
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


PARITY = {24:qr_parity,12:golay_parity,4:hamming_parity}


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
