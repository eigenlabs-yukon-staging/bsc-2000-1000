"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter, deque
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
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from dynamics import derive, data, fixed_tails
from qr import N, ORBITS, ROW_CAP, PERMUTATIONS, orbit_check, profiles as qr_compute

T=92
Q=1<<333
WEIGHTS=tuple(9**(665-t) for t in range(T+1))
DEN=10**665


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
    for weight in range(1,333):
        gap=333-2*weight
        prev,value=0,1
        for r,row in enumerate(signed):
            row.append(abs(value))
            if r<T:
                nxt,rem=divmod(gap*value-(332-r)*prev,r+1)
                assert rem==0
                prev,value=value,nxt
    result=prefix(665)[1:T+2]
    for k in range(1,GROUPS+1):
        capacities=[
            sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
            for w in range(1,333)
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
    return [2*d for d in result]


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


def weighted_miss(tails,alpha):
    anchor=[comb(M,w)*9**(M-w) if w<=M else 0 for w in range(T+1)]
    outside=[comb(M-1,w)*9**(M-1-w) if w<M else 0 for w in range(T+1)]
    remain=[comb(665-M,w)*9**(665-M-w) for w in range(T+1)]
    anchor_fail=Fraction(0)
    for d in range(RADIUS+1,min(M,T)+1):
        anchor_fail+=anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
    anchor_fail/=DEN
    diffs=[x-y for x,y in zip(alpha,alpha[1:]+[Fraction(0)],strict=True)]
    g=[]
    for s in range(T+1):
        p=sum(outside[w]*diffs[s+w] for w in range(T-s+1))/10**(M-1)
        g.append(upper_decimal(p,SCALE))
    ratios={}
    for den in NOISE_DENOMINATORS[1:]:
        ratios[den]=[(den**444*9**(444-z),10**444*(den-1)**(444-z))
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


RUN_GRID = (
    ((82,2),(83,4),(83,6),(83,7),(83,8),(84,8)),
    ((80,1),(80,2),(81,3),(81,4),(81,5),(81,6),(82,7),(82,8)),
)


def eighth_up(n,d,j):
    numerator=n**j<<512
    denominator=d**j
    floor=isqrt(isqrt(isqrt(numerator//denominator)))
    return floor+(floor**8*denominator<numerator)


@lru_cache(None)
def run_factors(first_punctured,a,j):
    qr,golay,repetition=[hist for _,hist in data()]
    shapes=[(48,qr,4),(24,golay,1),(2,repetition,2)]
    if first_punctured:
        shapes.append((1,{(1,1):1},1))
    else:
        shapes.append((2,repetition,1))
    factors=[]
    for bits,hist,count in shapes:
        row=[0]*(bits+1)
        denominator=100**(8*bits)
        for h,multiplicity in hist.items():
            numerator=sum(n*a**(8*w)*100**(8*(bits-w)) for w,n in enumerate(h))
            cap=eighth_up(numerator,denominator,j)
            for w,n in enumerate(h):
                row[w]+=multiplicity*n*cap
        factors.append((row,count))
    return factors


def run_upper(length,first_punctured,parity,free,a,j):
    paired=run_factors(first_punctured,a,j)
    coeff=product([*paired,(binomial(free),1)],T,4000)
    exponent,rem=divmod(j*free,8)
    top,bottom=100**8+a**8,100**8
    pow_two,frac_two=divmod(j*(parity-length),8)
    numerator=top**exponent*eighth_up(top,bottom,rem)*eighth_up(1,1<<frac_two,1)
    denominator=bottom**exponent*(1<<(pow_two+10*64))
    return [Fraction(n*numerator*100**(j*w),denominator*a**(j*w))
            for w,n in enumerate(coeff)]


def geometry(start,length):
    fields=[(start+i)%GROUPS for i in range(length+1)]
    first_punctured=start==0
    parity=length*M-(0 in fields[1:])
    free=(length-1)*M+parity
    local=free+2*M-first_punctured
    return first_punctured,parity,free,local


def loss(length,alpha,shapes):
    weighted=[w*a for w,a in zip(WEIGHTS,alpha,strict=True)]
    profiles=Counter(geometry(start,length) for start in range(GROUPS))
    result=Fraction(0)
    for (first_punctured,parity,free,local),multiplicity in profiles.items():
        bounds=[Fraction(comb(local,w)) for w in range(T+1)]
        for a,j in RUN_GRID[length-1]:
            estimate=run_upper(length,first_punctured,parity,free,a,j)
            bounds=[min(x,y) for x,y in zip(bounds,estimate,strict=True)]
        outside=binomial(665-local)
        result+=multiplicity*sum(b*sum(outside[t-w]*weighted[t] for t in range(w,T+1))
                                 for w,b in enumerate(bounds))/DEN
    return result


def certificate(data):
    previous = volume = success = 0
    alpha = []
    for t,bound in enumerate(full_direction_bound()):
        volume += comb(665,t)
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
        if i in (3,37):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert QR_ROW_BASE + ORBITS*ROW_CAP == BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert table_word(SQUARE_BASE)==0
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[12]+13)==PARITY[12](13)
    assert table_word(LIST_BASE[12]+4097)==candidate_row(12,1)[1]
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((12,0),(12,4095)):
        row = candidate_row(width,syndrome)
        parts = [(word>>17,(word>>12)&31,word&4095) for word in row]
        assert all(cost > 0 for cost, _, _ in parts)
        assert parts == sorted(parts)
        assert sorted(e for _,_,e in parts) == list(range(1<<width))
        for _,d,e in parts:
            assert d == e.bit_count()+(PARITY[width](e)^syndrome).bit_count()
    msg=(1<<1000)-37
    assert encode(msg)&((1<<1000)-1)==msg
    assert encode(msg).bit_length()==2000 and encode(msg)>>1996==15
    assert encode(0)==0


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
        w,_,nominal=derive(10)
        assert w == 3123857783472802875057418115050
        assert w**500 <= 1<<50651
        assert fixed_tails(nominal)==tails["10"]
        print(json.dumps({"W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
    elif part==2:
        for den in NOISE_DENOMINATORS[1:]:
            _,_,success=derive(den)
            assert fixed_tails(success)==tails[str(den)]
        print(json.dumps({"checked_tilts":list(NOISE_DENOMINATORS[1:])}))
    else:
        baseline,losses,alpha=certificate(shapes)
        anchor,child=weighted_miss(tails,alpha)
        component=baseline+sum(losses)+anchor+child
        assert 0<=component<=1
        extra_miss=Fraction(sum(comb(5,w)*9**(5-w) for w in range(3,6)),10**5)
        total=1-(1-extra_miss)*(1-component)**3
        assert total <= Fraction(487492753406772138,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "component_rounded":upper_decimal(component,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
