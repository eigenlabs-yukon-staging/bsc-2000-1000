"""Exact leader weights in one of four invariant syndrome sectors."""
from fractions import Fraction
from functools import reduce
from itertools import combinations
from operator import xor
from pathlib import Path
import json
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import BASE_SYNDROME_COLUMNS, TAIL_SYNDROME_COLUMNS, PARITY_ORDER, OFFSET, TAIL_LIMIT, permute_parity
from ops_samples import check_runtime



BASIS = (12656841, 1360143429, 337315020, 279144581, 379652293, 67314819, 3273785549, 3496149382, 107038536, 1353335370, 614632008, 11780932, 97616070, 1, 2, 4, 8, 16, 32, 64, 128, 256, 1024, 2048, 4096, 8192, 16384, 32768, 262144, 134217728, 2171086, 78782541)


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

    projections = [permute_parity(inverse(1 << i)) >> 24 for i in range(32)]

    def projections_for(bits):
        values = [0]
        for c in bits:
            values.extend(v ^ c for v in values[:])
        return values

    low_masks = [0] * 256
    for i, c in enumerate(projections_for(projections[:low])):
        low_masks[c] |= singleton[i]
    target = (OFFSET >> 24) ^ (permute_parity(inverse(sector << 30)) >> 24)
    forbidden = [low_masks[c ^ target] for c in
                 projections_for(projections[low:30])]
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
    count=excluded=0
    seen=[]
    blocked=[]
    for w in range(13):
        if w<=5:
            mark(rows,columns,w)
        else:
            limit=min(len(order),(33,36,40,44,len(order))[min(w-6,4)])
            new=rows.copy()
            expand(rows,new,order[:limit],permutations,low,full,count)
            mark(new,order[limit:]+remainder,w)
            rows=new
        n=sum(row.bit_count() for row in rows)
        bad=sum((row & mask).bit_count() for row, mask in zip(rows, forbidden))
        seen.append(n-count)
        blocked.append(bad-excluded)
        count,excluded=n,bad
        if count==size:
            break
    assert count==size
    return tuple(seen),tuple(blocked)



def tail_histogram():
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
    for w in range(25):
        if w <= 4:
            mark(rows, columns, w)
        else:
            limit = (32, 37, 42, 48, len(columns))[min(w - 5, 4)]
            new = rows.copy()
            expand(rows, new, columns[:limit], permutations, low, full, count)
            mark(new, columns[limit:], w)
            rows = new
        n = sum(row.bit_count() for row in rows)
        bad = sum((rows[row] & bits).bit_count() for row, bits in forbidden.items())
        seen.append(n - count)
        blocked.append(bad - excluded)
        count, excluded = n, bad
        if count == 1 << 24:
            break
    assert count == 1 << 24 and excluded == 48
    return tuple(seen), tuple(blocked)


WIDE_COUNTS = ((1, 46, 1080, 17466, 220416, 2312694, 20732246, 151973545, 605861899, 292434986, 187443, 2), (0, 6, 312, 8066, 138320, 1770606, 17965008, 142263338, 602483917, 308888783, 223468), (0, 6, 312, 8066, 138320, 1770606, 17963490, 142251435, 602677463, 308710973, 221153), (0, 6, 312, 8066, 138320, 1770606, 17965113, 142300977, 602937174, 308404437, 216813))
WIDE_EXCLUDED = ((0, 0, 0, 18, 486, 7038, 73109, 574348, 2371889, 1166680, 736, 0), (0, 0, 0, 7, 266, 5003, 61493, 532868, 2356658, 1237127, 882), (0, 0, 0, 4, 250, 5080, 62286, 534581, 2355617, 1235515, 971), (0, 0, 0, 9, 294, 5157, 61822, 533448, 2359226, 1233487, 861))
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
               for counts, excluded in zip(WIDE_COUNTS, WIDE_EXCLUDED))
    tail = safe_numerator(TAIL_COUNTS, TAIL_EXCLUDED, 56)
    p_wide = Fraction(wide * ((1 << 32) - TAIL_LIMIT), 10**64 * (1 << 32))
    p_tail = (Fraction(tail * ((1 << 32) - 48), 10**56 * (1 << 32))
              * Fraction(9, 10)**8)
    return p_wide**30 * p_tail


def main():
    sector = json.load(sys.stdin)["sector"]
    assert sector in (0, 1, 2, 3) and rank() == 32
    assert sorted(PARITY_ORDER) == list(range(32))
    hist, excluded = histogram(sector)
    assert hist == WIDE_COUNTS[sector]
    assert excluded == WIDE_EXCLUDED[sector]
    manifest = json.loads((HERE / "submission.json").read_text())
    error, exponent = [Fraction(int(c["score"]["numerator"]),
                                int(c["score"]["denominator"]))
                       for c in manifest["claims"]]
    assert 1 - combined_success() <= error
    data = {"sector": sector, "syndrome_rank": 32, "histogram": hist,
            "excluded": excluded, "B": str(error), "L": str(exponent),
            "error_passed": True}
    if sector == 0:
        tails, reserved = tail_histogram()
        assert tails == TAIL_COUNTS
        assert reserved == TAIL_EXCLUDED
        data.update(tail_histogram=tails, tail_excluded=reserved)
        data.update(check_runtime(exponent))
    print(json.dumps(data))


if __name__ == "__main__":
    main()
