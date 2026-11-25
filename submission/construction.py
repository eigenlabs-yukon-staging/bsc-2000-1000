"""Fixed field, parity cells and decoder RAM tables."""
from dataclasses import dataclass

from qr import N, ORBITS, PERMUTATIONS, ROW_CAP, orbit_labels, parity as qr_parity, permute_byte, row_words

M = 50
MAX_GROUPS = 3
WORDS = 1
HIGH_BITS = M - 64 * (WORDS - 1)
MASK = (1 << M) - 1
MASK64 = (1 << 64) - 1
BYTES = (M + 7) // 8
NIBBLES = (M + 3) // 4
MOD_LOW = 1 | (1 << 1) | (1 << 2) | (1 << 9)
MODULUS = (1 << M) | MOD_LOW
PRECISION = 16
LIMIT = 850
WIDTHS = (24,24,1,1)
OFFSETS = tuple(sum(WIDTHS[:j]) for j in range(len(WIDTHS)))
MAX_CELLS = (MAX_GROUPS - 1) * len(WIDTHS)
SQUARE_BYTES = tuple(
    sum(((value >> bit) & 1) << (2*bit) for bit in range(8)).to_bytes(2,"little")
    for value in range(256)
)


SAVED_B = 48
COMPONENT_GROUPS = (3, 3, 3, 3, 3, 3, 2)
FIELD_COUNT = sum(COMPONENT_GROUPS)
FIRST_WORK_RAM = SAVED_B + FIELD_COUNT*WORDS
ARRAYS = {name:FIRST_WORK_RAM+WORDS*i for i,name in enumerate(
    tuple(f"{label}{j}" for label in ("Y","P","X","C","B") for j in range(MAX_GROUPS))
    + ("SQ","AP"))}
NEXT = FIRST_WORK_RAM + WORDS*len(ARRAYS)
LAST = NEXT+M+1
SELECTED = LAST+M+1
FRAMES = SELECTED+MAX_CELLS
MULT = FRAMES+6*MAX_CELLS
REDUCE_1 = MULT+16*WORDS
REDUCE_4 = REDUCE_1+2
SQUARE_BASE = REDUCE_4+16
PARITY_BASE = {1:SQUARE_BASE+BYTES*256*WORDS}
LIST_BASE = {1:PARITY_BASE[1]+2}
QR_PARITY_BASE = LIST_BASE[1] + (1<<2)
QR_ORBIT_BASE = QR_PARITY_BASE + N
QR_PERM_BASE = QR_ORBIT_BASE + N
QR_ROW_BASE = QR_PERM_BASE + PERMUTATIONS * 6 * 256
QR_ROW_END = QR_ROW_BASE + ORBITS * ROW_CAP
PROFILE_COUNTS = {24:14,1:2}


@dataclass(frozen=True)
class Component:
    groups: int
    budgets: tuple
    cutoff: int
    budget_base: int
    profile24: int
    profile1: int
    cap_offsets: tuple

    @property
    def schedule(self):
        return WIDTHS*(self.groups-1)

    @property
    def radius(self):
        return len(self.budgets)-1

    @property
    def payload(self):
        return M*self.groups

    @property
    def length(self):
        return 2*self.payload


def component_layout(groups, budgets, cutoff, start):
    qr_profiles=start+len(budgets)
    unit_profiles=qr_profiles+ORBITS
    cap_base=unit_profiles+2
    caps=[cap_base]
    for width in WIDTHS*(groups-1):
        caps.append(caps[-1]+(LIMIT+1)*PROFILE_COUNTS[width])
    return Component(groups,budgets,cutoff,start,qr_profiles,unit_profiles,tuple(caps))


THREE=component_layout(3,(480, 433, 391, 356, 308, 269, 229, 188, 138, 92, 57, 8),48,QR_ROW_END)
TWO=component_layout(2,(473, 389, 372, 303, 264, 211, 177, 137, 76, 38),38,THREE.cap_offsets[-1])
MEMORY_END=TWO.cap_offsets[-1]
SPECS=(THREE,TWO)


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
        value=(value&MASK)^carryless(hi,MOD_LOW)
    return value


def square(value):
    octets=value.to_bytes((value.bit_length()+7)//8,"little")
    return remainder(int.from_bytes(b"".join(SQUARE_BYTES[i] for i in octets),"little"))


def cube(value):
    return remainder(carryless(value,square(value)))


PARITY = {24:qr_parity,1:int}


def linear(value):
    result = 0
    for w,i in zip(WIDTHS,OFFSETS,strict=True):
        result |= PARITY[w]((value>>i)&((1<<w)-1)) << i
    return result


def encode_component(message, groups):
    data=tuple((message>>(M*i))&MASK for i in range(groups))
    parity=tuple(linear(data[i])^cube(data[i-1]) for i in range(groups))
    return sum(value<<(M*i) for i,value in enumerate(data+parity))


def encode(message):
    result=0
    start=0
    for groups in COMPONENT_GROUPS:
        bits=M*groups
        low=(1<<bits)-1
        local=encode_component((message>>start)&low,groups)
        result|=(local&low)<<start
        result|=(local>>bits)<<(M*FIELD_COUNT+start)
        start+=bits
    return result


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
    if SQUARE_BASE<=address<PARITY_BASE[1]:
        row,word=divmod(address-SQUARE_BASE,WORDS)
        index,key=divmod(row,256)
        return (square_entry(index,key)>>(64*word))&MASK64
    if PARITY_BASE[1]<=address<LIST_BASE[1]:
        return address-PARITY_BASE[1]
    if LIST_BASE[1]<=address<QR_PARITY_BASE:
        syndrome,rank=divmod(address-LIST_BASE[1],2)
        return candidate_row(1,syndrome)[rank]
    if QR_PARITY_BASE<=address<QR_ORBIT_BASE:
        return qr_parity(address-QR_PARITY_BASE)
    if QR_ORBIT_BASE<=address<QR_PERM_BASE:
        return orbit_labels()[0][address-QR_ORBIT_BASE]
    if QR_PERM_BASE<=address<QR_ROW_BASE:
        row,value=divmod(address-QR_PERM_BASE,256)
        perm,byte=divmod(row,6)
        return permute_byte(perm,byte,value)
    if QR_ROW_BASE<=address<QR_ROW_END:
        orbit,ordinal=divmod(address-QR_ROW_BASE,ROW_CAP)
        return row_words(orbit,PRECISION)[ordinal]
    for spec in SPECS:
        if spec.budget_base<=address<spec.profile24:
            return spec.budgets[address-spec.budget_base]
        if spec.profile24<=address<spec.cap_offsets[-1]:
            from dynamics import adaptive_word
            return adaptive_word(spec,address)
    raise ValueError("not a fixed table address")
