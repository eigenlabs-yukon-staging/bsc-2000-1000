"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
from math import comb, isqrt
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (GROUPS,LIMIT,M,MEMORY_END,MODULUS,PARITY,PRECISION,
    RADIUS,REDUCE_4,BUDGETS,BUDGET_BASE,CAP_BASE,CAP_OFFSETS,PROFILE_BASE,TOTAL_CELLS,candidate_row,cost,encode,
    ANCHOR_BITS,QR_ROW_WORDS,QR_REP_BASE,QR_COST_BASE,poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE)

from dynamics import derive, data
from joint_admission import OUTPUT_SCALE, child_miss_upper
from qr import N, ORBITS, ROW_CAP, PERMUTATIONS, orbit_check, profiles as qr_compute

T=230
Q=1<<1000
WEIGHTS=tuple(9**(2000-t) for t in range(T+1))
DEN=10**2000


@lru_cache(None)
def binomial(n):
    return tuple(comb(n,i) for i in range(T+1))


@lru_cache(None)
def prefix(n):
    result=[0]
    for k in range(n+1):
        result.append(result[-1]+comb(n,k))
    return result


def ranked_groups(values):
    rows=sorted(values,reverse=True)
    counts=[0]
    magnitudes=[0]
    for value,count in rows:
        counts.append(counts[-1]+count)
        magnitudes.append(magnitudes[-1]+count*value)
    return rows,counts,magnitudes


def envelope(table,capacity):
    rows,counts,magnitudes=table
    i=bisect_left(counts,capacity)-1
    return magnitudes[i]+(capacity-counts[i])*rows[i][0]


def full_direction_bound():
    signed=[[] for _ in range(T+1)]
    for weight in range(1,1001):
        gap=1001-2*weight
        prev,value=0,1
        for r,row in enumerate(signed):
            row.append(abs(value))
            if r<T:
                nxt,rem=divmod(gap*value-(1000-r)*prev,r+1)
                assert rem==0
                prev,value=value,nxt
    result=prefix(2000)[1:T+2]
    for k in range(1,GROUPS+1):
        capacities=[
            3**k*sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
            for w in range(1,1001)
        ]
        assert all(c>=0 for c in capacities)
        sys=binomial(k*M)
        contributions=[0]*(T+1)
        for r,row in enumerate(signed):
            table=ranked_groups((v,c) for v,c in zip(row,capacities,strict=True) if c)
            for a in range(min(k*M,T-r)+1):
                contributions[a+r]+=envelope(table,sys[a])
        outside=binomial((GROUPS-k)*M)
        for t in range(T+1):
            result[t]+=comb(GROUPS,k)*sum(
                outside[b]*contributions[t-b] for b in range(min(t,(GROUPS-k)*M)+1)
            )
    assert all(a<=b for a,b in zip(result,result[1:]))
    return result


def product(parts, degree, bits):
    stride = (bits+7)//8
    mask = (1 << (8*stride*(degree+1))) - 1
    acc = 1
    for row, exponent in parts:
        encoded = int.from_bytes(b''.join(x.to_bytes(stride,'little') for x in row[:degree+1]),'little')
        while exponent:
            if exponent & 1:
                acc = (acc*encoded)&mask
            exponent >>= 1
            if exponent:
                encoded = (encoded*encoded)&mask
    buf = acc.to_bytes(stride*(degree+1),'little')
    return [int.from_bytes(buf[i:i+stride],'little') for i in range(0,len(buf),stride)]

def anchor_miss(alpha):
    anchor=[comb(ANCHOR_BITS,w)*9**(ANCHOR_BITS-w) if w<=ANCHOR_BITS else 0
            for w in range(T+1)]
    remain=[comb(2000-ANCHOR_BITS,w)*9**(2000-ANCHOR_BITS-w) for w in range(T+1)]
    miss=sum(anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
             for d in range(RADIUS+1,min(ANCHOR_BITS,T)+1))
    return miss/DEN


def root_up(n, d):
    scaled = n << 128
    s = isqrt(scaled//d)
    return s + (s*s*d < scaled)

def loss(length, a, b, alpha, data):
    free = M*(2*length-1)
    local = free+2*M
    parts = []
    for (width, hist), count in zip(data,(20,1,2),strict=True):
        row = [0]*(2*width+1)
        for profile, multiplicity in hist.items():
            h = sum(n*a**(2*i)*b**(4*width-2*i) for i,n in enumerate(profile))
            root = root_up(h,b**(4*width))
            for i,n in enumerate(profile):
                row[i] += multiplicity*n*root
        parts.append((row,count))
    parts.append(([comb(free,i) for i in range(T+1)],1))
    row = product(parts,T,5000)
    weighted = [n*b**i*a**(T-i) for i,n in enumerate(row)]
    outside = [comb(2000-local,i) for i in range(T+1)]
    shells = product(((weighted,1),(outside,1)),T,7000)
    scale = Fraction((a*a+b*b)**(free//2), b**free * (1<<(((M-1)*length)//2+23*64)) * a**T * DEN)
    return scale*sum(n*w*f for n,w,f in zip(shells,WEIGHTS,alpha))

def certificate(data):
    previous = volume = success = 0
    alpha = []
    for t,bound in enumerate(full_direction_bound()):
        volume += comb(2000,t)
        k = 1+bound//Q
        current = max(previous, volume*(2*k*Q-bound)//(Q*k*(k+1)))
        assert current<=Q
        success += (current-previous)*WEIGHTS[t]
        previous = current
        alpha.append(Fraction(2,k+1))
    ml0 = 1-Fraction(success,DEN)
    losses = [loss(length,a,100,alpha,data) for length,a in enumerate((42,),1)]
    return ml0, losses, alpha


def field_checks():
    x=2
    for i in range(1,M+1):
        x=square(x)
        if i in (100,250):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert QR_ROW_BASE+ORBITS*QR_ROW_WORDS==QR_REP_BASE
    assert QR_ROW_WORDS*2==ROW_CAP==N
    assert QR_REP_BASE+ORBITS==QR_COST_BASE
    assert QR_COST_BASE+14*49==BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((4,0),(4,15),(12,0),(12,2026),(12,4095)):
        row = candidate_row(width,syndrome)
        parts = [(word>>17,(word>>12)&31,word&4095) for word in row]
        assert all(cost > 0 for cost, _, _ in parts)
        assert parts == sorted(parts)
        assert sorted(e for _,_,e in parts) == list(range(1<<width))
        for _,d,e in parts:
            assert d == e.bit_count()+(PARITY[width](e)^syndrome).bit_count()
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg
    assert encode(msg).bit_length()<=2000


def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    shapes=data()
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values())}))
        return
    if part==1:
        field_checks()
        orbits,fixed=orbit_check()
        w,_,_=derive(10)
        assert w == 13247138559373527692065589430924434214824014876332140534580415970507848611677114108572430161516946984341966486032318854369464901
        assert w**1000 <= 1<<422291
        print(json.dumps({"W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
        return
    declared=json.loads((HERE/"joint.json").read_text())
    assert 0<=declared["child_miss_numerator"]<=OUTPUT_SCALE
    if part==2:
        alpha=[Fraction(2,2+d//Q) for d in full_direction_bound()]
        result=child_miss_upper(alpha)
        assert result==declared["child_miss_numerator"]
        print(json.dumps({"child_miss_numerator":result,"scale":OUTPUT_SCALE}))
    else:
        baseline,losses,alpha=certificate(shapes)
        anchor=anchor_miss(alpha)
        child=Fraction(declared["child_miss_numerator"],OUTPUT_SCALE)
        total=baseline+GROUPS*sum(losses)+anchor+child
        assert total <= Fraction(55175259435786241,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
