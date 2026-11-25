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


BASIS=(34049, 2097732, 142901504, 537976832, 1612709888, 86028, 1081369, 2117636, 4980738, 6307841, 17334912, 38010880, 134225929, 136314956, 268502112, 335544320, 538444160, 1074274384, 1091190784, 1241522304, 2164278272, 24, 518, 1072, 1408, 65929, 139268, 526400, 1048706, 1053056, 25821188, 143425536)


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
    assert len(order)==50 and len(remainder)==14

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


WIDE_COUNTS=((1, 50, 1251, 21000, 267381, 2762910, 23936952, 165804672, 605090562, 275595521, 261518, 6), (0, 4, 225, 6234, 113475, 1527120, 16119503, 131303370, 579662357, 344390461, 619067, 8), (0, 5, 270, 7215, 127260, 1667241, 17204953, 137280368, 587775829, 329192137, 486538, 8), (0, 5, 270, 7215, 127260, 1667241, 17204848, 137212238, 587204416, 329817894, 500429, 8))
EXCLUDED_COUNTS=((0, 0, 0, 0, 6, 127, 1346, 9139, 28400, 10151, 7, 0), (0, 0, 0, 0, 5, 100, 1031, 7334, 27260, 13428, 18, 0), (0, 0, 0, 1, 19, 192, 1430, 8551, 27759, 11191, 9, 0), (0, 0, 0, 2, 31, 250, 1544, 8392, 27153, 11767, 13, 0))

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
    assert sum(map(sum,EXCLUDED_COUNTS))==48+3*65536
    assert 1-Fraction(safe**31*tail,10**2000)<=B
    data={"sector":sector,"rank":32,"histogram":h,"excluded":e,
          "safe_main_numerator":safe,"tail_numerator":tail,"B":str(B),"L":str(L)}
    if sector==0:
        data.update(check_runtime(L))
    print(json.dumps(data))


if __name__=="__main__":
    main()
