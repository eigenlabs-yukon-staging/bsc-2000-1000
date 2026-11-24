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
    RADIUS,REDUCE_4,SLOPE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,
    ANCHOR_BITS,PAR_WIDTHS,SYS_WIDTHS,SEARCH_WIDTHS,SYS_ONLY_BASE)

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

def profiles(width, parity_map):
    weights = [i.bit_count() for i in range(1 << width)]
    values = [parity_map(i) for i in range(1 << width)]
    histogram = Counter()
    for syndrome in range(1 << width):
        profile = [0] * (2*width+1)
        for w, p in zip(weights, values):
            profile[w + weights[p ^ syndrome]] += 1
        histogram[tuple(profile)] += 1
    return histogram

DATA = [(w, profiles(w,PARITY[w])) for w in (12,4,1)]

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

def admission():
    parts=[]
    for (width,hist), exponent in zip(DATA,(54,4,2),strict=True):
        row=[0]*(LIMIT+1)
        weights=[9**(2*width-d) for d in range(2*width+1)]
        for profile,multiplicity in hist.items():
            normalizer=sum(a*b for a,b in zip(profile,weights,strict=True))
            for n,w in zip(profile,weights,strict=True):
                if n:
                    c=cost(normalizer,w)
                    if c<=LIMIT:
                        row[c]+=multiplicity*n*w
        parts.append((row,exponent))
    row=[0]*(LIMIT+1)
    for w in (9,1):
        row[cost(10,w)] += w
    parts.append((row,2))
    row=product(parts,LIMIT,4500)
    for i in range(1,len(row)):
        row[i]+=row[i-1]
    admitted=sum(comb(ANCHOR_BITS,i)*9**(ANCHOR_BITS-i)*row[LIMIT-SLOPE*i]
                 for i in range(RADIUS+1))
    return 1-Fraction(admitted,10**1666)


def root_up(n, d):
    scaled = n << 128
    s = isqrt(scaled//d)
    return s + (s*s*d < scaled)

SHAPES = (
    ((24,DATA[0][1],27),(8,DATA[1][1],2),(1,{(1,0):1,(0,1):1},2)),
    ((24,DATA[0][1],27),(8,DATA[1][1],2),(1,{(1,1):1},2)),
    ((24,DATA[0][1],27),(8,DATA[1][1],2),(2,DATA[2][1],2)),
)


@lru_cache(None)
def paired_parts(kind,a,b):
    parts=[]
    for bits,hist,count in SHAPES[kind]:
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
    scale_den=b**(2*half_free)*(1<<((parity_bits-length-2*includes_anchor)//2+31*64))
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
    assert len(SEARCH_WIDTHS)==TOTAL_CELLS==62
    assert MEMORY_END<=1<<32
    for e in (0,1):
        assert table_word(SYS_ONLY_BASE+e) == e|(e<<12)|(cost(10,9**(1-e))<<17)
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg and encode(msg).bit_length()<=2000


def operation_bound():
    profile_cost={}
    for width,hist in DATA:
        rows=[]
        for profile in hist:
            masses=[9**(2*width-d) for d in range(2*width+1)]
            normalizer=sum(n*w for n,w in zip(profile,masses))
            counts=Counter()
            for n,w in zip(profile,masses):
                if n:
                    counts[cost(normalizer,w)]+=n
            rows.append(tuple(sorted(counts.items())))
        profile_cost[width]=rows
    profile_cost[0]=[tuple((cost(10,9**(1-e)),1) for e in (0,1))]
    k=WORDS
    shift=3*k+3
    square=BYTES*(5+3*k)+k
    multiply=3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops=square+multiply
    anchor = 410
    leaf = 250
    widths=SEARCH_WIDTHS
    dp=[0]*(LIMIT+1)
    for i in range(len(widths)-1,-1,-1):
        width=widths[i]
        reserve=len(widths)-1-i
        inc=28
        if i==len(widths)-1:inc+=cube_ops+leaf
        elif (i+1)%len(WIDTHS)==0:inc+=cube_ops+180
        out=[0]*(LIMIT+1)
        incrow=[k+inc for k in dp]
        for row in profile_cost[width]:
            w=[0]*(LIMIT+1)
            for cc,n in row:
                for rem in range(cc+reserve,LIMIT+1):
                    w[rem]+=n*incrow[rem-cc]
            out=[max(x,y) for x,y in zip(out,w)]
        dp=[40+x for x in out]
    return 1000+sum(comb(ANCHOR_BITS,dep)*(44+anchor+cube_ops+dp[LIMIT-SLOPE*dep]) for dep in range(RADIUS+1))

def main():
    part=json.load(sys.stdin)["part"]
    assert part in range(2)
    runs=((), (1, 2))
    result={}
    if runs[part]:
        baseline,alpha=certificate()
        base=upper_decimal(baseline,10**18)
        assert base <= 55181057551535181
        losses=[upper_decimal(length_charge(a,alpha),10**18) for a in runs[part]]
        assert all(v <= (373970, 395514)[a-1] for a,v in zip(runs[part],losses,strict=True))
        result.update({"baseline":base,"losses":losses})
    if part==0:
        field_checks()
        miss=upper_decimal(admission(),10**18)
        assert miss <= 1997643287089
        w=operation_bound()
        assert w == 1714363491736672601245918726983446372521051114591435656093119217826676119531910802664750335410883572974736894
        assert w**500 <= 1<<179773
        result.update({"miss":miss,"W":str(w),"memory_end":MEMORY_END})
    print(json.dumps(result))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
