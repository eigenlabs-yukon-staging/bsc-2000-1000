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
from construction import SYNDROME_COLUMNS, TAIL_SYNDROME_COLUMNS, OFFSET, TAIL_LIMIT
from ops_samples import check_runtime



BASIS=(30722, 172490752, 303300736, 1, 5144, 33792, 266496, 1050630, 21234176, 268828704, 285999106, 287834176, 570441731, 1073973248, 1078988800, 1275142144, 2214600896, 2232418304, 2684452866, 16, 2440, 6656, 16400, 38912, 74000, 98944, 131264, 135444, 401408, 8388866, 33564688, 1077970960)


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
    columns=tuple(transform(c) for c in SYNDROME_COLUMNS)
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

    projections = [inverse(1 << i) >> 24 for i in range(32)]

    def projections_for(bits):
        values = [0]
        for c in bits:
            values.extend(v ^ c for v in values[:])
        return values

    low_masks = [0] * 256
    for i, c in enumerate(projections_for(projections[:low])):
        low_masks[c] |= singleton[i]
    target = (OFFSET >> 24) ^ (inverse(sector << 30) >> 24)
    forbidden = [low_masks[c ^ target] for c in
                 projections_for(projections[low:30])]
    order=tuple(sorted((c for c in columns if c>>30==0),
                       key=lambda v:(v&mask).bit_count()))
    remainder=tuple(c for c in columns if c>>30!=0)
    assert len(order)==49 and len(remainder)==15

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


WIDE_COUNTS = ((1, 49, 1207, 20063, 254533, 2635637, 22991559, 161401337, 604411546, 281783431, 242455, 6), (0, 5, 269, 7171, 126323, 1654393, 17079642, 136528126, 587849558, 330045239, 451090, 8), (0, 6, 314, 8152, 140108, 1794514, 18163664, 142385497, 594795766, 316093062, 360733, 8), (0, 4, 226, 6278, 114412, 1539968, 16249176, 132361159, 583555273, 339410434, 504886, 8))
WIDE_EXCLUDED = ((0, 0, 1, 42, 802, 9888, 91373, 644609, 2370235, 1076492, 862, 0), (0, 0, 1, 29, 510, 6673, 68517, 542531, 2305396, 1269050, 1597, 0), (0, 0, 0, 16, 449, 6701, 70917, 560469, 2330063, 1224312, 1377, 0), (0, 0, 0, 8, 289, 5136, 60764, 515690, 2292294, 1318237, 1886, 0))
TAIL_COUNTS = (1, 56, 1540, 27710, 363368, 3350166, 11197863, 1836511, 1)
TAIL_EXCLUDED = (0, 0, 0, 0, 1, 9, 33, 5, 0)


def safe_numerator(weight_counts, excluded, length):
    assert len(weight_counts) == len(excluded)
    assert all(0 <= b <= n for n, b in zip(weight_counts, excluded))
    return sum((n - b) * 9**(length - w)
               for w, (n, b) in enumerate(zip(weight_counts, excluded)))


def hamming_lane():
    columns = (3, 5, 6, 7, 1, 2, 4)
    success = 0
    for error in range(128):
        bits = [(error >> i) & 1 for i in range(7)]
        a, b, c, d, p0, p1, p2 = bits
        l0, l1, l2 = a ^ d, b ^ d, c ^ d
        s0, s1, s2 = l0 ^ b ^ p0, l0 ^ c ^ p1, l1 ^ c ^ p2
        s = s0 | (s1 << 1) | (s2 << 2)
        e = reduce(xor, (v for v, bit in zip(columns, bits) if bit), 0)
        assert s == e
        output = (l0 ^ (s0 & s1), l1 ^ (s0 & s2), l2 ^ (s1 & s2), d ^ (s0 & s1 & s2))
        assert (not any(output)) == (error.bit_count() <= 1)
        if not any(output):
            success += 9**(7 - error.bit_count())
    return Fraction(success, 10**7)


def combined_success(hamming):
    assert sum(map(sum, WIDE_COUNTS)) == 1 << 32
    assert sum(map(sum, WIDE_EXCLUDED)) == TAIL_LIMIT
    assert sum(TAIL_COUNTS) == TAIL_LIMIT and sum(TAIL_EXCLUDED) == 48
    wide = sum(safe_numerator(counts, excluded, 64)
               for counts, excluded in zip(WIDE_COUNTS, WIDE_EXCLUDED))
    tail = safe_numerator(TAIL_COUNTS, TAIL_EXCLUDED, 56)
    p_wide = Fraction(wide * ((1 << 32) - TAIL_LIMIT), 10**64 * (1 << 32))
    p_tail = (Fraction(tail * ((1 << 32) - 48), 10**56 * (1 << 32))
              * Fraction(9, 10)**8)
    return hamming**192 * p_wide**6 * p_tail


def main():
    sector = json.load(sys.stdin)["sector"]
    assert sector in (0, 1, 2, 3) and rank() == 32
    hist, excluded = histogram(sector)
    assert hist == WIDE_COUNTS[sector]
    assert excluded == WIDE_EXCLUDED[sector]
    manifest = json.loads((HERE / "submission.json").read_text())
    error, exponent = [Fraction(int(c["score"]["numerator"]),
                                int(c["score"]["denominator"]))
                       for c in manifest["claims"]]
    hamming = hamming_lane()
    assert hamming == Fraction(531441, 625000)
    assert 1 - combined_success(hamming) <= error
    data = {"sector": sector, "syndrome_rank": 32, "histogram": hist,
            "excluded": excluded, "B": str(error), "L": str(exponent),
            "error_passed": True}
    if sector == 0:
        tails, reserved = tail_histogram()
        assert tails == TAIL_COUNTS
        assert reserved == TAIL_EXCLUDED
        data.update(tail_histogram=tails, tail_excluded=reserved,
                    hamming_success=str(hamming))
        data.update(check_runtime(exponent))
    print(json.dumps(data))


if __name__ == "__main__":
    main()
