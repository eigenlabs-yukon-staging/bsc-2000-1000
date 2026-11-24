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
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

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

def data():
    qr_data=json.loads((HERE/"profiles.json").read_text())
    hist={tuple(row["h"]):row["multiplicity"] for row in qr_data}
    assert len(hist)==14 and sum(hist.values())==N
    return [(24,hist)]+[(w,profiles(w,PARITY[w])) for w in (4,1)]


def row_counts(width,profile):
    cap=ROW_CAP if width==24 else 1<<width
    result=[]
    for n in profile:
        included=min(n,cap)
        result.append(included)
        cap-=included
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

SCALE=10**60
NOISE_DENOMINATORS=(10,12)


def child_tails(width_profiles,den):
    parts=[]
    for (width,hist),exponent in zip(width_profiles,(30,6,6),strict=True):
        row=[0]*(LIMIT+1)
        nominal=[9**(2*width-d) for d in range(2*width+1)]
        tilted=[(den-1)**(2*width-d) for d in range(2*width+1)]
        for profile,multiplicity in hist.items():
            normalizer=sum(a*b for a,b in zip(profile,nominal,strict=True))
            for n,w,q in zip(row_counts(width,profile),nominal,tilted,strict=True):
                if n:
                    c=cost(normalizer,w)
                    if c<=LIMIT:
                        row[c]+=multiplicity*n*q
        parts.append((row,exponent))
    norm=den**1500
    row=product(parts,LIMIT,norm.bit_length())
    for i in range(1,len(row)):
        row[i]+=row[i-1]
    tails=[]
    for d in range(RADIUS+1):
        miss=norm-row[LIMIT-SLOPE*d]
        assert 0<=miss<=norm
        tails.append((SCALE*miss+norm-1)//norm)
    return tails


def weighted_miss(tails,alpha):
    anchor=[comb(M,w)*9**(M-w) if w<=M else 0 for w in range(T+1)]
    outside=anchor
    remain=[comb(2000-M,w)*9**(2000-M-w) for w in range(T+1)]
    anchor_fail=Fraction(0)
    for d in range(RADIUS+1,min(M,T)+1):
        anchor_fail+=anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
    anchor_fail/=DEN
    diffs=[x-y for x,y in zip(alpha,alpha[1:]+[Fraction(0)],strict=True)]
    g=[]
    for s in range(T+1):
        p=sum(outside[w]*diffs[s+w] for w in range(T-s+1))/10**M
        g.append(upper_decimal(p,SCALE))
    ratios={}
    for den in NOISE_DENOMINATORS[1:]:
        ratios[den]=[(den**1500*9**(1500-z),10**1500*(den-1)**(1500-z))
                     for z in range(T+1)]
    child_fail=0
    for d in range(RADIUS+1):
        subtotal=0
        for z in range(T-d+1):
            best=tails["10"][d]
            for den in NOISE_DENOMINATORS[1:]:
                numerator,denom=ratios[den][z]
                adjusted=(tails[str(den)][d]*numerator+denom-1)//denom
                best=min(best,adjusted)
            subtotal+=g[d+z]*best
        child_fail+=anchor[d]*subtotal
    return anchor_fail,Fraction(child_fail,10**M*SCALE**2)


def root_up(n, d):
    scaled = n << 128
    s = isqrt(scaled//d)
    return s + (s*s*d < scaled)

def run_coefficients(length, a, b, data):
    free = M*(2*length-1)
    parts = []
    for (width, hist), count in zip(data,(10,2,2),strict=True):
        row = [0]*(2*width+1)
        for profile, multiplicity in hist.items():
            h = sum(n*a**(2*i)*b**(4*width-2*i) for i,n in enumerate(profile))
            root = root_up(h,b**(4*width))
            for i,n in enumerate(profile):
                row[i] += multiplicity*n*root
        parts.append((row,count))
    parts.append(([comb(free,i) for i in range(T+1)],1))
    row = product(parts,T,4000)
    scale_num = (a*a+b*b)**(free//2)
    scale_den = b**free * (1<<(((M-1)*length)//2+14*64))
    return [Fraction(n*scale_num*b**w,scale_den*a**w) for w,n in enumerate(row)]


def loss(length, alpha, data):
    local = (2*length+1)*M
    by_weight = [Fraction(comb(local,w)) for w in range(T+1)]
    for a in range(28,55,2):
        by_weight = [min(x,y) for x,y in zip(by_weight,run_coefficients(length,a,100,data),strict=True)]
    outside = binomial(2000-local)
    tilted = [w*g for w,g in zip(WEIGHTS,alpha,strict=True)]
    propagate = [sum(outside[t-w]*tilted[t] for t in range(w,T+1)) for w in range(T+1)]
    return sum(x*y for x,y in zip(by_weight,propagate,strict=True))/DEN


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
    losses = [loss(length,alpha,data) for length in range(1,GROUPS)]
    return ml0, losses, alpha


def field_checks():
    x=2
    for i in range(1,M+1):
        x=square(x)
        if i in (50,125):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert QR_ROW_BASE + ORBITS*ROW_CAP == MEMORY_END<=1<<32
    assert table_word(SQUARE_BASE)==0
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[4]+13)==PARITY[4](13)
    assert table_word(LIST_BASE[4]+21)==candidate_row(4,1)[5]
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((4,0),(4,15)):
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


def operation_bound(data):
    profile_cost={}
    for width,hist in data:
        rows=[]
        for profile in hist:
            masses=[9**(2*width-d) for d in range(2*width+1)]
            normalizer=sum(n*w for n,w in zip(profile,masses,strict=True))
            counts=Counter()
            for n,w in zip(row_counts(width,profile),masses,strict=True):
                if n:
                    counts[cost(normalizer,w)]+=n
            rows.append(tuple(sorted(counts.items())))
        profile_cost[width]=rows
    k=WORDS
    shift=3*k+3
    square=BYTES*(5+3*k)+k
    multiply=3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops=square+multiply
    anchor = 250
    leaf = 200
    widths=WIDTHS*(GROUPS-1)
    dp=[0]*(LIMIT+1)
    for i in range(len(widths)-1,-1,-1):
        width=widths[i]
        reserve=len(widths)-1-i
        inc=80 if width==24 else 28
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
        overhead=60 if width==24 else 40
        dp=[overhead+x for x in out]
    return 1000+sum(comb(M,dep)*(44+anchor+cube_ops+dp[LIMIT-SLOPE*dep]) for dep in range(RADIUS+1))

def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    shapes=data()
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values())}))
        return
    tails=json.loads((HERE/"tails.json").read_text())
    assert set(tails)=={str(q) for q in NOISE_DENOMINATORS}
    assert all(len(row)==RADIUS+1 and all(0<=x<=SCALE for x in row) for row in tails.values())
    if part==1:
        field_checks()
        orbits,fixed=orbit_check()
        assert child_tails(shapes,10)==tails["10"]
        w=operation_bound(shapes)
        assert w == 17453369834132354728616933631134082069953296005161442271916533470499301182969893511286
        assert w**125 <= 1<<35396
        print(json.dumps({"W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
    elif part==2:
        for den in NOISE_DENOMINATORS[1:]:
            assert child_tails(shapes,den)==tails[str(den)]
        print(json.dumps({"checked_tilts":list(NOISE_DENOMINATORS[1:])}))
    else:
        baseline,losses,alpha=certificate(shapes)
        anchor,child=weighted_miss(tails,alpha)
        total=baseline+GROUPS*sum(losses)+anchor+child
        assert total <= Fraction(55205024520759224,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
