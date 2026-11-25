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
from construction import SYNDROME_COLUMNS, OFFSET, TAIL_BASE, tail_leaders
from ops_samples import check_runtime


BASIS=(30722, 172490752, 303300736, 1, 5144, 33792, 266496, 1050630, 21234176, 268828704, 285999106, 287834176, 570441731, 1073973248, 1078988800, 1275142144, 2214600896, 2232418304, 2684452866, 16, 2440, 6656, 16400, 38912, 74000, 98944, 131264, 135444, 401408, 8388866, 33564688, 1077970960)


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
    for addresses in (range(48),range(TAIL_BASE,1<<32)):
        for address in addresses:
            s=address^OFFSET
            z=tables[0][s&255]^tables[1][s>>8&255]^tables[2][s>>16&255]^tables[3][s>>24]
            if z>>30==sector:
                row=(z&local_mask)>>low
                forbidden[row]=forbidden.get(row,0)|singleton[z&mask]
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
        bad=sum((rows[j]&v).bit_count() for j,v in forbidden.items())
        seen.append(n-count)
        blocked.append(bad-excluded)
        count,excluded=n,bad
        if count==size:
            break
    assert count==size
    return tuple(seen),tuple(blocked)



WIDE_COUNTS=((1, 49, 1207, 20063, 254533, 2635637, 22991559, 161401337, 604411546, 281783431, 242455, 6), (0, 5, 269, 7171, 126323, 1654393, 17079642, 136528126, 587849558, 330045239, 451090, 8), (0, 6, 314, 8152, 140108, 1794514, 18163664, 142385497, 594795766, 316093062, 360733, 8), (0, 4, 226, 6278, 114412, 1539968, 16249176, 132361159, 583555273, 339410434, 504886, 8))
EXCLUDED_COUNTS=((0, 0, 0, 0, 0, 6, 133, 1588, 8812, 5839, 6, 0), (0, 0, 0, 0, 0, 3, 104, 1458, 8734, 6099, 2, 0), (0, 0, 0, 0, 0, 1, 76, 1299, 8639, 6387, 14, 0), (0, 0, 0, 0, 0, 1, 60, 1074, 8051, 7192, 6, 0))

def main():
    sector=json.load(sys.stdin)["sector"]
    assert sector in (0,1,2,3) and rank()==32
    h,e=histogram(sector)
    assert h==WIDE_COUNTS[sector] and e==EXCLUDED_COUNTS[sector]
    d=json.loads((HERE/"submission.json").read_text())
    B,L=[Fraction(int(x["score"]["numerator"]),int(x["score"]["denominator"]))
         for x in d["claims"]]
    safe=sum((n-bad)*9**(64-w)
             for h0,e0 in zip(WIDE_COUNTS,EXCLUDED_COUNTS)
             for w,(n,bad) in enumerate(zip(h0,e0)))
    tail=sum(9**(16-error.bit_count()) for error in tail_leaders())
    assert sum(map(sum,EXCLUDED_COUNTS))==48+65536
    data_safe=(1<<32)-48-(1<<16)
    assert 1-Fraction((safe*data_safe)**31*tail,10**2000*(1<<992))<=B
    data={"sector":sector,"rank":32,"histogram":h,"excluded":e,
          "safe_main_numerator":safe,"tail_numerator":tail,"B":str(B),"L":str(L)}
    if sector==0:
        data.update(check_runtime(L))
    print(json.dumps(data))


if __name__=="__main__":
    main()
