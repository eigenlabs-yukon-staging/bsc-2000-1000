"""Exact leader weights in one of four invariant syndrome sectors."""
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations
from operator import xor
from pathlib import Path
import json
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import BASE_SYNDROME_COLUMNS, TAIL_SYNDROME_COLUMNS, PARITY_ORDER, OFFSET, TAIL_LIMIT, MAIN_COLUMNS, RESTORE_MAX_WIDE, permute_parity
from ops_samples import check_runtime



BASIS = (12656841, 1360143429, 337315020, 279144581, 379652293, 67314819, 3273785549, 3496149382, 107038536, 1353335370, 614632008, 11780932, 97616070, 1, 2, 8, 16, 32, 128, 256, 4096, 134217728, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 67108864, 2171086, 78782541)


def expand(rows, new, selected, permutations, low, full, count):
    mask = (1 << low) - 1
    edges = [(c >> low, permutations[c & mask]) for c in selected]
    if 4 * count >= 3 * (len(rows) << low):
        for j, union in enumerate(rows):
            if union == full:
                continue
            for high, swaps in edges:
                row = rows[j ^ high]
                for shift, left, right in swaps:
                    row = ((row & left) << shift) | ((row & right) >> shift)
                union |= row
                if union == full:
                    break
            new[j] = union
    else:
        for high, swaps in edges:
            if not swaps:
                for j, row in enumerate(rows):
                    new[j ^ high] |= row
                continue
            for j, row in enumerate(rows):
                if row:
                    for shift, left, right in swaps:
                        row = ((row & left) << shift) | ((row & right) >> shift)
                    new[j ^ high] |= row


def rank():
    pivots={}
    for item in BASIS:
        assert 0<item<1<<32
        a=item
        while a:
            i=a.bit_length()-1
            if i not in pivots:
                pivots[i]=a
                break
            a^=pivots[i]
    return len(pivots)


def transform(v):
    return sum(((v&row).bit_count()&1)<<j for j,row in enumerate(BASIS))


def histogram(sector):
    columns=tuple(transform(c) for c in BASE_SYNDROME_COLUMNS)
    low=13
    mask=(1<<low)-1
    full=(1<<(1<<low))-1
    swap_masks=[]
    for j in range(low):
        shift=1<<j
        left=full//((1<<shift)+1)
        swap_masks.append((shift,left,left<<shift))
    permutations=[tuple(swap_masks[j] for j in range(low) if v>>j&1)
                  for v in range(1<<low)]
    singleton=tuple(1<<j for j in range(1<<low))
    size = 1 << 30
    local_mask = size - 1
    # Invert the map by eliminating pairs (transformed, original).
    pivots = {}
    for i in range(32):
        orig = 1 << i
        vec = transform(orig)
        while vec:
            j = vec.bit_length() - 1
            if j not in pivots:
                pivots[j] = (vec, orig)
                break
            a, b = pivots[j]
            vec ^= a
            orig ^= b

    def inverse(z):
        v = 0
        while z:
            a, b = pivots[z.bit_length() - 1]
            z ^= a
            v ^= b
        return v

    assert BASIS[22:30] == tuple(1 << i for i in PARITY_ORDER[24:])
    weights = tail_statistics()[2]
    tail_types = len(RESTORE_MAX_WIDE)

    low_values = [0]
    for i in range(low):
        delta = permute_parity(inverse(1 << i))
        assert delta < TAIL_LIMIT
        low_values.extend(a ^ delta for a in low_values[:])

    parity_bytes = []
    for byte in range(3):
        values = [0]
        for i in range(8):
            delta = MAIN_COLUMNS[8 * byte + i] & 255
            values.extend(a ^ delta for a in values[:])
        parity_bytes.append(values)

    marked_rows = []
    for middle in range(1 << 9):
        index = ((OFFSET >> 24) << 9) | middle
        syndrome = permute_parity(inverse((sector << 30) | (index << low)))
        assert syndrome >> 24 == OFFSET >> 24
        base = syndrome ^ OFFSET
        groups = [0] * tail_types
        zeros = [0] * tail_types
        for i, delta in enumerate(low_values):
            key = delta ^ base
            if key < 48:
                continue
            tail_weight = weights[key ^ (OFFSET >> 8)]
            groups[tail_weight] |= singleton[i]
            check = (parity_bytes[0][key & 255] ^
                     parity_bytes[1][key >> 8 & 255] ^
                     parity_bytes[2][key >> 16])
            if check == OFFSET & 255:
                zeros[tail_weight] |= singleton[i]
        marked_rows.append((index, groups, zeros))
    order=tuple(sorted((c for c in columns if c>>30==0),
                       key=lambda v:(v&mask).bit_count()))
    remainder=tuple(c for c in columns if c>>30!=0)
    assert len(order)+len(remainder)==64

    def mark(rows,choices,w):
        for group in combinations(choices,w):
            syndrome=reduce(xor,group,0)
            if syndrome>>30==sector:
                rows[(syndrome&local_mask)>>low]|=singleton[syndrome&mask]

    rows=[0]*(size>>low)
    count=0
    seen=[]
    blocked=[]
    joint=[]
    zero_check=[]
    for w in range(13):
        before = [rows[index] for index, _, _ in marked_rows]
        if w<=5:
            mark(rows,columns,w)
        else:
            limit=min(len(order),(33,36,40,44,len(order))[min(w-6,4)])
            new=rows.copy()
            expand(rows,new,order[:limit],permutations,low,full,count)
            mark(new,order[limit:]+remainder,w)
            rows=new
        n=sum(row.bit_count() for row in rows)
        low_count = 0
        joint_at_weight = [0] * tail_types
        zero_at_weight = [0] * tail_types
        for previous, (index, groups, zeros) in zip(before, marked_rows, strict=True):
            added = rows[index] ^ previous
            low_count += added.bit_count()
            for i, (group, zero_group) in enumerate(zip(groups, zeros, strict=True)):
                joint_at_weight[i] += (added & group).bit_count()
                zero_at_weight[i] += (added & zero_group).bit_count()
        seen.append(n-count)
        blocked.append(low_count)
        joint.append(tuple(joint_at_weight))
        zero_check.append(tuple(zero_at_weight))
        count=n
        if count==size:
            break
    assert count==size
    return tuple(seen),tuple(blocked),tuple(joint),tuple(zero_check)



@lru_cache(None)
def tail_statistics():
    low = 10
    mask = (1 << low) - 1
    full = (1 << (1 << low)) - 1
    swap_masks = []
    for j in range(low):
        shift = 1 << j
        left = full // ((1 << shift) + 1)
        swap_masks.append((shift, left, left << shift))
    permutations = [tuple(swap_masks[j] for j in range(low) if v >> j & 1)
                    for v in range(1 << low)]
    singleton = tuple(1 << j for j in range(1 << low))
    columns = tuple(sorted(TAIL_SYNDROME_COLUMNS,
                           key=lambda c: (c & mask).bit_count()))
    rows = [0] * (1 << (24 - low))
    forbidden = {}
    for r in range(48):
        z = r ^ (OFFSET >> 8)
        index = z >> low
        forbidden[index] = forbidden.get(index, 0) | singleton[z & mask]

    def mark(rows, choices, weight):
        for locations in combinations(choices, weight):
            s = reduce(xor, locations, 0)
            rows[s >> low] |= singleton[s & mask]

    count = excluded = 0
    seen = []
    blocked = []
    weights = bytearray(1 << 24)
    for w in range(25):
        before = rows.copy()
        if w <= 4:
            mark(rows, columns, w)
        else:
            limit = (32, 37, 42, 48, len(columns))[min(w - 5, 4)]
            new = rows.copy()
            expand(rows, new, columns[:limit], permutations, low, full, count)
            mark(new, columns[limit:], w)
            rows = new
        for index, (old, row) in enumerate(zip(before, rows, strict=True)):
            change = row ^ old
            while change:
                one = change & -change
                weights[(index << low) | (one.bit_length() - 1)] = w
                change ^= one
        n = sum(row.bit_count() for row in rows)
        bad = sum((rows[row] & bits).bit_count() for row, bits in forbidden.items())
        seen.append(n - count)
        blocked.append(bad - excluded)
        count, excluded = n, bad
        if count == 1 << 24:
            break
    assert count == 1 << 24 and excluded == 48
    return tuple(seen), tuple(blocked), weights


WIDE_COUNTS = ((1, 46, 1080, 17466, 220416, 2312694, 20732246, 151973545, 605861899, 292434986, 187443, 2), (0, 6, 312, 8066, 138320, 1770606, 17965008, 142263338, 602483917, 308888783, 223468), (0, 6, 312, 8066, 138320, 1770606, 17963490, 142251435, 602677463, 308710973, 221153), (0, 6, 312, 8066, 138320, 1770606, 17965113, 142300977, 602937174, 308404437, 216813))
WIDE_EXCLUDED = ((0, 0, 0, 18, 486, 7038, 73109, 574348, 2371889, 1166680, 736, 0), (0, 0, 0, 7, 266, 5003, 61493, 532868, 2356658, 1237127, 882), (0, 0, 0, 4, 250, 5080, 62286, 534581, 2355617, 1235515, 971), (0, 0, 0, 9, 294, 5157, 61822, 533448, 2359226, 1233487, 861))
WIDE_JOINT = (((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 4, 12, 2), (0, 0, 0, 0, 5, 98, 333, 50), (0, 0, 1, 13, 160, 1343, 4848, 673), (0, 1, 9, 139, 1572, 15153, 49329, 6906), (0, 2, 62, 965, 12898, 116926, 390854, 52639), (0, 8, 221, 3990, 51275, 487967, 1610540, 217864), (0, 3, 111, 1919, 25562, 238366, 793519, 107178), (0, 0, 0, 0, 24, 137, 512, 63), (0, 0, 0, 0, 0, 0, 0, 0)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 2, 5, 0), (0, 0, 0, 0, 8, 57, 176, 25), (0, 0, 1, 7, 101, 1031, 3398, 465), (0, 0, 6, 108, 1322, 12520, 41706, 5831), (0, 0, 51, 903, 11552, 108288, 363064, 49010), (0, 8, 227, 3877, 51096, 484102, 1601202, 216146), (1, 14, 131, 2059, 27108, 253717, 842098, 111999), (0, 0, 0, 0, 13, 179, 605, 85)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 2, 1, 1), (0, 0, 0, 0, 9, 46, 169, 26), (0, 0, 0, 13, 90, 1018, 3522, 437), (0, 0, 4, 78, 1384, 12747, 42186, 5887), (0, 2, 41, 827, 11696, 108993, 364438, 48584), (0, 7, 195, 4002, 50666, 485369, 1600079, 215299), (0, 3, 136, 1963, 27461, 252721, 840440, 112791), (0, 0, 0, 1, 22, 196, 663, 89)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 2, 7, 0), (0, 0, 0, 0, 7, 63, 190, 34), (0, 0, 0, 12, 93, 1010, 3581, 461), (0, 1, 9, 94, 1333, 12765, 41869, 5751), (0, 2, 59, 839, 11788, 109622, 362501, 48637), (0, 2, 172, 3974, 51183, 485847, 1601834, 216214), (0, 3, 104, 1936, 27374, 252168, 838973, 112929), (0, 0, 0, 1, 19, 193, 561, 87)))
WIDE_ZERO_CHECK = (((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 1), (0, 0, 0, 0, 1, 4, 24, 0), (0, 0, 0, 1, 0, 51, 195, 26), (0, 0, 0, 5, 54, 435, 1574, 223), (0, 0, 0, 15, 202, 2020, 6232, 841), (0, 0, 0, 8, 109, 923, 3025, 413), (0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 0, 0, 0)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 4, 12, 2), (0, 0, 0, 1, 3, 50, 180, 28), (0, 0, 0, 4, 37, 425, 1409, 184), (0, 0, 0, 20, 203, 1845, 6219, 875), (0, 0, 0, 12, 111, 989, 3310, 455), (0, 0, 0, 0, 0, 1, 4, 0)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 5, 11, 2), (0, 0, 0, 0, 10, 50, 178, 28), (0, 0, 0, 5, 45, 424, 1386, 184), (0, 0, 0, 16, 214, 1883, 6194, 831), (0, 0, 0, 13, 107, 1003, 3349, 440), (0, 0, 0, 0, 0, 1, 4, 1)), ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 2, 17, 0), (0, 0, 0, 0, 7, 45, 149, 21), (0, 0, 0, 3, 43, 420, 1429, 188), (0, 0, 1, 14, 207, 1871, 6278, 854), (0, 0, 0, 13, 114, 987, 3298, 416), (0, 0, 0, 0, 0, 0, 5, 1)))
TAIL_COUNTS = (1, 56, 1540, 27720, 365821, 3442653, 11403249, 1536176)
TAIL_EXCLUDED = (0, 0, 0, 0, 0, 1, 34, 13)


def safe_numerator(weight_counts, excluded, length):
    assert len(weight_counts) == len(excluded)
    assert all(0 <= b <= n for n, b in zip(weight_counts, excluded))
    return sum((n - b) * 9**(length - w)
               for w, (n, b) in enumerate(zip(weight_counts, excluded)))


def combined_success():
    assert sum(map(sum, WIDE_COUNTS)) == 1 << 32
    assert sum(map(sum, WIDE_EXCLUDED)) == TAIL_LIMIT
    assert sum(TAIL_COUNTS) == TAIL_LIMIT and sum(TAIL_EXCLUDED) == 48
    wide = sum(safe_numerator(counts, excluded, 64)
               for counts, excluded in zip(WIDE_COUNTS, WIDE_EXCLUDED, strict=True))
    tail = safe_numerator(TAIL_COUNTS, TAIL_EXCLUDED, 56)
    first_keys = (1 << 32) - TAIL_LIMIT
    tail_keys = [0] * len(TAIL_COUNTS)
    for sector, zeros in zip(WIDE_JOINT, WIDE_ZERO_CHECK, strict=True):
        for wide_weight, (classes, zero_classes) in enumerate(zip(sector, zeros, strict=True)):
            for tail_weight, (count, zero_count) in enumerate(zip(classes, zero_classes, strict=True)):
                assert 0 <= zero_count <= count
                tail_keys[tail_weight] += count
                if wide_weight <= RESTORE_MAX_WIDE[tail_weight]:
                    first_keys += count
                    wide += count * 9**(64 - wide_weight)
                    tail -= count * 9**(56 - tail_weight)
                else:
                    first_keys += zero_count
    assert tail_keys == [total - bad for total, bad in zip(TAIL_COUNTS, TAIL_EXCLUDED, strict=True)]
    p_wide = Fraction(wide * first_keys, 10**64 * (1 << 32))
    p_tail = (Fraction(tail * ((1 << 32) - 48), 10**56 * (1 << 32))
              * Fraction(9, 10)**8)
    return p_wide**30 * p_tail


def main():
    sector = json.load(sys.stdin)["sector"]
    assert sector in (0, 1, 2, 3) and rank() == 32
    assert sorted(PARITY_ORDER) == list(range(32))
    hist, excluded, joint, zero_check = histogram(sector)
    assert hist == WIDE_COUNTS[sector]
    assert excluded == WIDE_EXCLUDED[sector]
    assert joint == WIDE_JOINT[sector]
    assert zero_check == WIDE_ZERO_CHECK[sector]
    manifest = json.loads((HERE / "submission.json").read_text())
    error, exponent = [Fraction(int(c["score"]["numerator"]),
                                int(c["score"]["denominator"]))
                       for c in manifest["claims"]]
    assert 1 - combined_success() <= error
    data = {"sector": sector, "syndrome_rank": 32, "histogram": hist,
            "excluded": excluded, "joint": joint, "zero_check": zero_check,
            "B": str(error), "L": str(exponent),
            "error_passed": True}
    tails, reserved, _ = tail_statistics()
    assert tails == TAIL_COUNTS and reserved == TAIL_EXCLUDED
    if sector == 0:
        data.update(tail_histogram=tails, tail_excluded=reserved)
        data.update(check_runtime(exponent))
    print(json.dumps(data))


if __name__ == "__main__":
    main()
