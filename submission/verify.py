"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
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

T=59
Q=1<<200
WEIGHTS=tuple(9**(400-t) for t in range(T+1))
DEN=10**400


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
    for weight in range(1,201):
        gap=201-2*weight
        prev,value=0,1
        for r,row in enumerate(signed):
            row.append(abs(value))
            if r<T:
                nxt,rem=divmod(gap*value-(200-r)*prev,r+1)
                assert rem==0
                prev,value=value,nxt
    result=prefix(400)[1:T+2]
    for k in range(1,GROUPS+1):
        capacities=[
            3**k*sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
            for w in range(1,201)
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


SCALE=10**60
NOISE_DENOMINATORS=(10,12)


def weighted_miss(tails,alpha):
    anchor=[comb(M,w)*9**(M-w) if w<=M else 0 for w in range(T+1)]
    outside=anchor
    remain=[comb(400-M,w)*9**(400-M-w) for w in range(T+1)]
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
        ratios[den]=[(den**300*9**(300-z),10**300*(den-1)**(300-z))
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


@lru_cache(None)
def shifted_max(row,m):
    q=deque()
    out=[]
    # row index w-1 represents mask weight w. For sum S of ordinary masks,
    # maximize over first weight 1..m, which gives row[S:S+m].
    for i,value in enumerate(row):
        while q and q[-1][1]<=value:q.pop()
        q.append((i,value))
        if i>=m-1:
            while q[0][0]<=i-m:q.popleft()
            out.append(q[0][1])
    return out

@lru_cache(None)
def kernel(length):
    parity=M*length
    ordinary=M*(length-1)
    ordinary_mult=3 if M%2==0 else 1
    volume=prefix(parity+ordinary)
    result=[volume[min(t+1,len(volume)-1)]+volume[min(t,len(volume)-1)] for t in range(T+1)]
    rows=[[] for _ in range(T+1)]
    for weight in range(1,parity+1):
        gap=parity+1-2*weight
        previous,current=0,1
        for radius,row in enumerate(rows):
            row.append(abs(previous+current))
            if radius<T:
                following,rem=divmod(gap*current-(parity-radius)*previous,radius+1)
                assert rem==0
                previous,current=current,following
    for selected in range(length):
        inside_bits=selected*M
        outside_bits=ordinary-inside_bits
        capacities=[ordinary_mult**selected*sum(
            (-1)**(selected-j)*comb(selected,j)*comb(j*M,w)
            for j in range(1,selected+1)
        ) for w in range(1,inside_bits+1)]
        coefficients=[0]*(T+1)
        for radius,row in enumerate(rows):
            if not selected:
                coefficients[radius]+=max(row[:M])
                continue
            shifted=shifted_max(tuple(row),M)
            excluded=ranked_groups((row[w-1],count) for w,count in enumerate(capacities,1) if count)
            included=ranked_groups((shifted[w],count) for w,count in enumerate(capacities,1) if count)
            for w in range(min(inside_bits,T-radius)+1):
                cap=comb(inside_bits,w)
                coefficients[radius+w]+=envelope(excluded,cap)+envelope(included,cap)
        choose=comb(length-1,selected)
        for t in range(T+1):
            result[t]+=choose*sum(comb(outside_bits,b)*coefficients[t-b]
                                  for b in range(min(t,outside_bits)+1))
    return result


def conv(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                out[i+j]+=x*y
    return tuple(out)


@lru_cache(None)
def pair_profiles():
    result=Counter({(1,):1})
    shapes=dict(data())
    for width in WIDTHS:
        out=Counter()
        for left,count_left in result.items():
            for right,count_right in shapes[width].items():
                out[conv(left,right)]+=count_left*count_right
        result=out
    return result


def loss(length,alpha,data):
    parity_bits=M*length
    free_bits=M*(2*length-1)
    outside_bits=2*GROUPS*M-2*M-free_bits
    out=binomial(outside_bits)
    free=binomial(free_bits)
    direction=kernel(length)
    denom=1<<(parity_bits+1)
    bad=[0]*(T+1)
    for profile,multiplicity in pair_profiles().items():
        for weight in range(T+1):
            upper=min(weight,len(profile)-1)
            total=sum(profile[d]*free[weight-d] for d in range(upper+1))
            possible=sum(profile[d]*direction[weight-d] for d in range(upper+1))
            same=sum(profile[d]*free[weight-d]*direction[weight-d] for d in range(upper+1))
            bad[weight]+=multiplicity*min(denom*total,possible*total-same)
    weight_alpha=[a*w for a,w in zip(alpha,WEIGHTS,strict=True)]
    return sum(bad[w]*sum(out[t-w]*weight_alpha[t] for t in range(w,T+1))
               for w in range(T+1))/(denom*DEN)



def certificate(data):
    previous = volume = success = 0
    alpha = []
    for t,bound in enumerate(full_direction_bound()):
        volume += comb(400,t)
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
        if i in (10,25):
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
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((1,0),(1,1)):
        row = candidate_row(width,syndrome)
        parts = [(word>>17,(word>>12)&31,word&4095) for word in row]
        assert all(cost > 0 for cost, _, _ in parts)
        assert parts == sorted(parts)
        assert sorted(e for _,_,e in parts) == list(range(1<<width))
        for _,d,e in parts:
            assert d == e.bit_count()+(PARITY[width](e)^syndrome).bit_count()
    msg=(1<<1000)-37
    assert encode(msg)&((1<<1000)-1)==msg
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
    tails=json.loads((HERE/"tails.json").read_text())
    assert set(tails)=={str(q) for q in NOISE_DENOMINATORS}
    assert all(len(row)==RADIUS+1 and all(0<=x<=SCALE for x in row) for row in tails.values())
    if part==1:
        field_checks()
        orbits,fixed=orbit_check()
        w,_,nominal=derive(10)
        assert w == 2520497841559917180
        assert w**1000 <= 1<<61129
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
        component=baseline+GROUPS*sum(losses)+anchor+child
        assert 0<=component<=1
        total=1-(1-component)**5
        assert total <= Fraction(773356200462581483,10**18)
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
