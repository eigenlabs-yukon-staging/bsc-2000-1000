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
from construction import SYNDROME_COLUMNS, OFFSET, LOW_BASE, tail_leaders
from ops_samples import check_runtime

BASIS=(2147484242, 67240010, 1073872928, 1075839504, 2149581632, 45070, 65683, 2097265, 8519680, 21495840, 33557042, 33754112, 38043648, 1, 2, 4, 8, 16, 32, 64, 128, 512, 2048, 4096, 8192, 262144, 524288, 134217728, 268435456, 536870912, 1057825, 17413)
WIDE_COUNTS=((1, 51, 1298, 22073, 283220, 2932655, 25172450, 163491450, 467663405, 370424359, 43750535, 324, 3), (0, 5, 270, 7215, 127260, 1667241, 17128648, 130651175, 453122709, 404544238, 66492021, 1038, 4), (0, 5, 270, 7215, 127260, 1667241, 17129667, 130703846, 453381624, 404491612, 66232087, 993, 4), (0, 3, 178, 5161, 97636, 1357375, 14676076, 118003826, 442989786, 417502528, 79107232, 2019, 4))
EXCLUDED_COUNTS=((0, 0, 0, 1, 19, 189, 1385, 7871, 21312, 16521, 1866, 0, 0), (0, 0, 0, 1, 16, 144, 1069, 6724, 20783, 17713, 2714, 0, 0), (0, 0, 0, 0, 4, 86, 953, 6735, 20933, 17761, 2692, 0, 0), (0, 0, 0, 0, 4, 75, 794, 5887, 20459, 18620, 3325, 0, 0))

def expand(rows, new, selected, permutations, low, full, count):
    mask = (1 << low) - 1
    edges = [(c >> low, permutations[c & mask]) for c in selected]
    if 4 * count >= 3 * (1 << 30):
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
    unit=tuple(transform(1<<i) for i in range(32))
    tables=[]
    for byte in range(4):
        vals=[0]
        for j in range(8):
            vals.extend(x^unit[8*byte+j] for x in vals[:])
        tables.append(vals)
    size=1<<30
    local_mask=size-1
    forbidden={}
    for addresses in (range(48),range(LOW_BASE,1<<32)):
        for address in addresses:
            s=address^OFFSET
            z=tables[0][s&255]^tables[1][s>>8&255]^tables[2][s>>16&255]^tables[3][s>>24]
            if z>>30==sector:
                row=(z&local_mask)>>low
                forbidden[row]=forbidden.get(row,0)|singleton[z&mask]
    order=tuple(sorted((c for c in columns if c>>30==0),
                       key=lambda v:(v&mask).bit_count()))
    remainder=tuple(c for c in columns if c>>30!=0)
    assert len(order)==51 and len(remainder)==13

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
            limit=(33,36,40,44,51)[min(w-6,4)]
            new=rows.copy()
            expand(rows,new,order[:limit],permutations,low,full,count)
            mark(new,order[limit:]+remainder,w)
            rows=new
        n=sum(row.bit_count() for row in rows)
        bad=sum((rows[j]&v).bit_count() for j,v in forbidden.items())
        seen.append(n-count)
        blocked.append(bad-excluded)
        count,excluded=n,bad
    assert count==size
    return tuple(seen),tuple(blocked)


def main():
    sector=json.load(sys.stdin)["sector"]
    assert rank()==32
    h,e=histogram(sector)
    assert h==WIDE_COUNTS[sector] and e==EXCLUDED_COUNTS[sector]
    d=json.loads((HERE/"submission.json").read_text())
    B,L=[Fraction(int(x["score"]["numerator"]),int(x["score"]["denominator"]))
         for x in d["claims"]]
    safe=sum((n-bad)*9**(64-w)
             for h0,e0 in zip(WIDE_COUNTS,EXCLUDED_COUNTS)
             for w,(n,bad) in enumerate(zip(h0,e0)))
    tail=sum(9**(16-error.bit_count()) for error in tail_leaders())
    assert sum(map(sum,EXCLUDED_COUNTS))==48+3*65536
    assert 1-Fraction(safe**31*tail,10**2000)<=B
    data={"sector":sector,"rank":32,"histogram":h,"excluded":e,
          "safe_main_numerator":safe,"tail_numerator":tail,"B":str(B),"L":str(L)}
    if sector==0:
        data.update(check_runtime(L))
    print(json.dumps(data))


if __name__=="__main__":
    main()
