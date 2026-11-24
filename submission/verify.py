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
    RADIUS,REDUCE_4,ANCHOR_CAPS,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS)

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

DATA = [(w, profiles(w,PARITY[w])) for w in (12,4)]

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

def admission(precision, limit, radius):
    parts = []
    for (width, hist), exponent in zip(DATA,(41,2)):
        row = [0]*(limit+1)
        weights = [9**(2*width-d) for d in range(2*width+1)]
        for profile, multiplicity in hist.items():
            normalizer = sum(a*b for a,b in zip(profile,weights))
            for n, weight in zip(profile,weights):
                if n:
                    c = cost(normalizer,weight)
                    assert precision == PRECISION and c > 0
                    if c <= limit:
                        row[c] += multiplicity*n*weight
        parts.append((row,exponent))
    row = product(parts,limit,3400)
    for i in range(1,len(row)):
        row[i] += row[i-1]
    admitted = sum(comb(M,i)*9**(M-i)*row[ANCHOR_CAPS[i]] for i in range(radius+1))
    return 1-Fraction(admitted,10**1500)

def root_up(n, d):
    scaled = n << 128
    s = isqrt(scaled//d)
    return s + (s*s*d < scaled)

def loss(length, a, b, alpha):
    free = M*(2*length-1)
    local = free+2*M
    parts = []
    for (width, hist), count in zip(DATA,(41,2)):
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
    scale = Fraction((a*a+b*b)**(free//2), b**free * (1<<(((M-1)*length)//2+43*64)) * a**T * DEN)
    return scale*sum(n*w*f for n,w,f in zip(shells,WEIGHTS,alpha))

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
    losses = [loss(length,a,100,alpha) for length,a in enumerate((42,),1)]
    return ml0, losses, ml0+GROUPS*sum(losses)


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
    assert MEMORY_END<=1<<32
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
    k=WORDS
    shift=3*k+3
    square=BYTES*(5+3*k)+k
    multiply=3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops=square+multiply
    anchor = 560
    leaf = 320
    widths=WIDTHS*(GROUPS-1)
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
    return 1000+sum(comb(M,dep)*(44+anchor+cube_ops+dp[ANCHOR_CAPS[dep]]) for dep in range(RADIUS+1))

def main():
    mode=json.load(sys.stdin)["part"]
    assert mode in (0,1)
    if mode==0:
        field_checks()
        error=admission(PRECISION,LIMIT,RADIUS)
        assert error <= Fraction(1998443260,10**18)
        w = operation_bound()
        assert w == 4392293791981942423011995978321459076161109293240437652321034300500504129044423156112161897943266241236676774849669917403555714562176279168
        assert w**500 <= 1<<230281
        print(json.dumps({"admission_rounded":upper_decimal(error,10**18),
            "W":str(w),"memory_end":MEMORY_END}))
    else:
        baseline, losses, ml = certificate()
        assert ml <= Fraction(55175259434787520,10**18)
        print(json.dumps({"full_rounded":upper_decimal(ml,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses]}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
