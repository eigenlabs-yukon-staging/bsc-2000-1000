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
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,
    ANCHOR_BITS,PAR_WIDTHS,SYS_WIDTHS,SEARCH_WIDTHS,SYS_ONLY_BASE,
    QR_ROW_BASE,SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

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
    masks=[[] for _ in range(T+1)]
    for w in range(1,1001):
        prev,cur=0,1
        for r,row in enumerate(masks):
            row.append(abs(cur))
            if r<T:
                nxt,rem=divmod(cur*(1001-2*w)-(1000-r)*prev,r+1)
                assert rem==0
                prev,cur=cur,nxt
    result=prefix(2000)[1:T+2]
    for special in (False,True):
        for k in range(0 if special else 1,GROUPS):
            if special:
                cap=[16*3**k*sum((-1)**(k-h)*comb(k,h)*(comb(M*h+ANCHOR_BITS,w)-comb(M*h,w))
                            for h in range(k+1)) for w in range(1,1001)]
            else:
                cap=[3**k*sum((-1)**(k-h)*comb(k,h)*comb(M*h,w)
                            for h in range(k+1)) for w in range(1,1001)]
            n=M*k+ANCHOR_BITS*special
            inside=binomial(n)
            outside=binomial(1000-n)
            con=[0]*(T+1)
            for r,row in enumerate(masks):
                grouped=ranked_groups((v,c) for v,c in zip(row,cap,strict=True) if c)
                for a in range(min(n,T-r)+1):
                    con[a+r]+=envelope(grouped,(7 if special else 1)*inside[a])
            for t in range(T+1):
                result[t]+=comb(GROUPS-1,k)*sum(outside[b]*con[t-b] for b in range(min(t,1000-n)+1))
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

@lru_cache(None)
def paired_parts(kind,a,b):
    [(w0,qr),(w1,golay),(w2,hamming),(w3,repetition)]=data()
    common=((48,qr,13),(24,golay,1),(8,hamming,2))
    extras=(((1,{(1,0):1,(0,1):1},2),),((1,{(1,1):1},2),),((2,repetition,2),))
    shape=common+extras[kind]
    parts=[]
    for bits,hist,count in shape:
        row=[0]*(bits+1)
        for profile,multiplicity in hist.items():
            h=sum(n*a**(2*i)*b**(2*bits-2*i) for i,n in enumerate(profile))
            root=root_up(h,b**(2*bits))
            for i,n in enumerate(profile):
                row[i]+=multiplicity*n*root
        parts.append((row,count))
    return product(parts,T,4000)


def run_coefficients(length, kind, parity_bits, free, includes_anchor, a, b):
    paired=paired_parts(kind,a,b)
    row=product(((paired,1),(binomial(free),1)),T,4000)
    half_free=(free+1)//2
    scale_num=(a*a+b*b)**half_free
    scale_den=b**(2*half_free)*(1<<((parity_bits-length-2*includes_anchor)//2+18*64))
    return [Fraction(n*scale_num*b**w,scale_den*a**w) for w,n in enumerate(row)]


def geometry(start,length):
    fields=[(start+i)%GROUPS for i in range(length+1)]
    parity=sum(PAR_WIDTHS[i] for i in fields[1:])
    free=parity+sum(SYS_WIDTHS[i] for i in fields[1:-1])
    local=free+SYS_WIDTHS[start]+PAR_WIDTHS[start]
    return min(start,2),parity,free,local,0 in fields[:-1]


def length_charge(length,alpha):
    terms=Counter(geometry(i,length) for i in range(GROUPS))
    tilted=[w*g for w,g in zip(WEIGHTS,alpha,strict=True)]
    result=Fraction(0)
    for (kind,parity,free,local,includes_anchor), mult in terms.items():
        by_weight=[Fraction(comb(local,w)) for w in range(T+1)]
        for a in range(28,55,2):
            tilt=run_coefficients(length,kind,parity,free,includes_anchor,a,100)
            by_weight=[min(x,y) for x,y in zip(by_weight,tilt,strict=True)]
        outside=binomial(2000-local)
        propagate=[sum(outside[t-w]*tilted[t] for t in range(w,T+1)) for w in range(T+1)]
        result+=mult*sum(x*y for x,y in zip(by_weight,propagate,strict=True))/DEN
    return result


def certificate():
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
    return ml0, alpha

def field_checks():
    x=2
    for j in range(1,M+1):
        x=square(x)
        if j in (2,167):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert sum(SYS_WIDTHS)==sum(PAR_WIDTHS)==1000
    assert len(SEARCH_WIDTHS)==TOTAL_CELLS==36
    assert QR_ROW_BASE+ORBITS*ROW_CAP==BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[12]+1023)==PARITY[12](1023)
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    for e in (0,1):
        assert table_word(SYS_ONLY_BASE+e) == e|(e<<12)|(cost(10,9**(1-e))<<17)
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg and encode(msg).bit_length()<=2000


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
        assert w == 3160110756777792829726456158289815335633316720255983875772330598236424830974228163195405920592926
        assert w**500 <= 1<<160283
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
        baseline,alpha=certificate()
        losses=[length_charge(a,alpha) for a in range(1,GROUPS)]
        anchor=anchor_miss(alpha)
        child=Fraction(declared["child_miss_numerator"],OUTPUT_SCALE)
        total=baseline+sum(losses)+anchor+child
        assert total <= Fraction(55181127530032615,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
