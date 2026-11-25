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
    ANCHOR_BITS,PAR_WIDTHS,SYS_WIDTHS,SEARCH_WIDTHS,SYS_ONLY_BASE,
    QR_ROW_BASE,SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from kernel import kernel
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
                cap=[8*sum((-1)**(k-h)*comb(k,h)*(comb(M*h+ANCHOR_BITS,w)-comb(M*h,w))
                            for h in range(k+1)) for w in range(1,1001)]
            else:
                cap=[sum((-1)**(k-h)*comb(k,h)*comb(M*h,w)
                            for h in range(k+1)) for w in range(1,1001)]
            n=M*k+ANCHOR_BITS*special
            inside=binomial(n)
            outside=binomial(1000-n)
            con=[0]*(T+1)
            for r,row in enumerate(masks):
                grouped=ranked_groups((v,c) for v,c in zip(row,cap,strict=True) if c)
                for a in range(min(n,T-r)+1):
                    con[a+r]+=envelope(grouped,(3 if special else 1)*inside[a])
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

@lru_cache(None)
def data():
    qr_data=json.loads((HERE/"profiles.json").read_text())
    hist={tuple(row["h"]):row["multiplicity"] for row in qr_data}
    assert len(hist)==14 and sum(hist.values())==N
    return [(24,hist)]+[(w,profiles(w,PARITY[w])) for w in (12,4,1)]


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
    for (width,hist),exponent in zip(width_profiles,(30,6,12,17),strict=True):
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
    row=[0]*(LIMIT+1)
    row[cost(10,9)]=den-1
    row[cost(10,1)]=1
    parts.append((row,1))
    norm=den**1715
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
    outside=[comb(M,w)*9**(M-w) if w<=M else 0 for w in range(T+1)]
    anchor=[comb(ANCHOR_BITS,w)*9**(ANCHOR_BITS-w) if w<=ANCHOR_BITS else 0 for w in range(T+1)]
    remain=[comb(2000-ANCHOR_BITS,w)*9**(2000-ANCHOR_BITS-w) for w in range(T+1)]
    anchor_fail=Fraction(0)
    for d in range(RADIUS+1,min(ANCHOR_BITS,T)+1):
        anchor_fail+=anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
    anchor_fail/=DEN
    diffs=[x-y for x,y in zip(alpha,alpha[1:]+[Fraction(0)],strict=True)]
    g=[]
    for s in range(T+1):
        p=sum(outside[w]*diffs[s+w] for w in range(T-s+1))/10**M
        g.append(upper_decimal(p,SCALE))
    ratios={}
    for den in NOISE_DENOMINATORS[1:]:
        ratios[den]=[(den**1715*9**(1715-z),10**1715*(den-1)**(1715-z))
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
    return anchor_fail,Fraction(child_fail,10**ANCHOR_BITS*SCALE**2)


GRID_ROWS=((1, 0, 142, 142, 427, True, ((83, 4), (83, 6), (82, 1), (83, 8), (83, 1), (83, 7))), (1, 1, 143, 143, 428, False, ((83, 8), (83, 5), (83, 3), (83, 6))), (1, 2, 143, 143, 429, False, ((83, 3), (83, 5), (83, 8), (83, 7))), (2, 0, 285, 428, 713, True, ((81, 4), (80, 1), (81, 5), (81, 3), (80, 2), (81, 1), (81, 6))), (2, 1, 286, 429, 714, False, ((80, 2), (81, 4), (80, 1), (81, 3), (81, 5), (81, 6))), (2, 2, 285, 427, 713, True, ((80, 2), (81, 4), (80, 1), (81, 3), (81, 5), (81, 1), (81, 6))), (2, 2, 286, 429, 715, False, ((81, 3), (81, 4), (80, 1), (80, 2), (81, 5), (81, 1), (81, 6))), (3, 0, 428, 714, 999, True, ((79, 1), (79, 2), (80, 3), (80, 4), (80, 1), (80, 5))), (3, 1, 429, 715, 1000, False, ((79, 2), (79, 1), (80, 3), (80, 4), (80, 1), (80, 5))), (3, 2, 428, 713, 999, True, ((79, 2), (79, 1), (80, 3), (80, 4), (80, 1), (80, 5))), (3, 2, 429, 715, 1001, False, ((79, 2), (79, 1), (80, 3), (80, 4), (80, 1), (80, 5))), (4, 0, 571, 1000, 1285, True, ((79, 2), (79, 3), (79, 1), (79, 4), (79, 5))), (4, 1, 572, 1001, 1286, False, ((79, 3), (79, 1), (79, 2), (79, 4), (79, 5))), (4, 2, 571, 999, 1285, True, ((79, 3), (79, 1), (79, 2), (79, 4), (79, 5))), (4, 2, 572, 1001, 1287, False, ((79, 1), (79, 3), (79, 2), (79, 4), (79, 5))), (5, 0, 714, 1286, 1571, True, ((79, 3), (79, 2), (78, 1), (79, 4), (79, 5))), (5, 1, 715, 1287, 1572, False, ((79, 3), (79, 2), (79, 4), (78, 2), (78, 1), (79, 5))), (5, 2, 714, 1285, 1571, True, ((79, 3), (79, 2), (79, 4), (79, 1), (78, 1), (79, 5))), (5, 2, 715, 1287, 1573, False, ((79, 3), (79, 2), (79, 4), (78, 2), (78, 1), (79, 1), (79, 5))), (6, 0, 857, 1572, 1857, True, ((78, 4), (78, 3), (79, 5), (78, 5), (79, 6))), (6, 1, 858, 1573, 1858, False, ((78, 4), (78, 3), (79, 5), (78, 5), (79, 6))), (6, 2, 857, 1571, 1857, True, ((78, 4), (79, 5), (78, 3), (79, 4), (78, 5), (79, 6))))
GRID={(length,kind,p,u,n,special):params for length,kind,p,u,n,special,params in GRID_ROWS}


def power_up(n,d,j):
    numerator=n**j<<512
    denominator=d**j
    initial=isqrt(isqrt(isqrt(numerator//denominator)))
    return initial+(initial**8*denominator<numerator)


@lru_cache(None)
def paired_parts(kind,a,b,j):
    (_,qr),(_,golay),(_,hamming),(_,repetition)=data()
    common=((48,qr,5),(24,golay,1),(8,hamming,2))
    extras=(((2,repetition,2),(1,{(1,0):1,(0,1):1},1)),
            ((2,repetition,2),(1,{(1,1):1},1)),
            ((2,repetition,3),))
    parts=[]
    for bits,hist,count in common+extras[kind]:
        row=[0]*(bits+1)
        d=b**(8*bits)
        for profile,multiplicity in hist.items():
            n=sum(v*a**(8*i)*b**(8*bits-8*i) for i,v in enumerate(profile))
            root=power_up(n,d,j)
            for i,v in enumerate(profile):
                row[i]+=multiplicity*v*root
        parts.append((row,count))
    return product(parts,T,4000)


def run_coefficients(length,kind,parity_bits,free,includes_anchor,a,b,j):
    paired=paired_parts(kind,a,b,j)
    row=product(((paired,1),(binomial(free),1)),T,4000)
    exponent,rem=divmod(j*free,8)
    top,bottom=a**8+b**8,b**8
    e,rem2=divmod(j*(parity_bits-length-includes_anchor),8)
    numerator=top**exponent*power_up(top,bottom,rem)*power_up(1,1<<rem2,1)
    denominator=bottom**exponent*(1<<(e+(11+2)*64))
    plain=Fraction(numerator,denominator)
    projected_scale,projected_rem=divmod(j*parity_bits,8)
    scale=Fraction(power_up(1,1<<projected_rem,1),1<<(projected_scale+(11+2)*64))
    cumulative=Fraction(0)
    projected=Fraction(0)
    out=[]
    direction=projected_kernel(length,kind==0,includes_anchor and kind!=0)
    for w,n in enumerate(row):
        candidate=Fraction(direction[w]*a**(8*w),b**(8*w))
        if candidate>cumulative:
            cumulative=candidate
            projected=scale*power_up(cumulative.numerator,cumulative.denominator,j)
        out.append(min(plain,projected)*Fraction(n*b**(j*w),a**(j*w)))
    return out


@lru_cache(None)
def projected_kernel(length,first_short,other_short):
    return kernel(M,1,length,first_short,other_short)


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
        for a,j in GRID[length,kind,parity,free,local,includes_anchor]:
            tilt=run_coefficients(length,kind,parity,free,includes_anchor,a,100,j)
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
        if j in (11,13):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert sum(SYS_WIDTHS)==sum(PAR_WIDTHS)==1000
    assert len(SEARCH_WIDTHS)==TOTAL_CELLS==66
    assert QR_ROW_BASE+ORBITS*ROW_CAP==MEMORY_END<=1<<32
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[12]+1023)==PARITY[12](1023)
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    for e in (0,1):
        assert table_word(SYS_ONLY_BASE+e) == e|(e<<12)|(cost(10,9**(1-e))<<17)
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg and encode(msg).bit_length()<=2000


def operation_bound(data):
    profile_cost={}
    for width,hist in data:
        rows=[]
        for profile in hist:
            masses=[9**(2*width-d) for d in range(2*width+1)]
            normalizer=sum(n*w for n,w in zip(profile,masses))
            counts=Counter()
            for n,w in zip(row_counts(width,profile),masses,strict=True):
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
    anchor = 250
    leaf = 200
    widths=SEARCH_WIDTHS
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
    return 1000+sum(comb(ANCHOR_BITS,dep)*(44+anchor+cube_ops+dp[LIMIT-SLOPE*dep]) for dep in range(RADIUS+1))

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
        assert w == 14287706924187778761802264616515855042789125897848060586917352
        assert w**1000 <= 1<<203153
        print(json.dumps({"W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
    elif part==2:
        for den in NOISE_DENOMINATORS[1:]:
            assert child_tails(shapes,den)==tails[str(den)]
        print(json.dumps({"checked_tilts":list(NOISE_DENOMINATORS[1:])}))
    else:
        baseline,alpha=certificate()
        losses=[length_charge(a,alpha) for a in range(1,GROUPS)]
        anchor,child=weighted_miss(tails,alpha)
        total=baseline+sum(losses)+anchor+child
        assert total <= Fraction(68051660172330392,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
