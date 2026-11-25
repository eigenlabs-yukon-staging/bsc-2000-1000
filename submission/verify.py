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
    ANCHOR_BITS,COMPONENTS,RADIUS,REDUCE_4,BUDGETS,BUDGET_BASE,CAP_BASE,CAP_OFFSETS,PROFILE_BASE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from dynamics import derive, data
from joint_admission import OUTPUT_SCALE, child_miss_upper
from qr import N, ORBITS, ROW_CAP, PERMUTATIONS, orbit_check, profiles as qr_compute

T=38
Q=1<<100
WEIGHTS=tuple(9**(200-t) for t in range(T+1))
DEN=10**200


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
    for weight in range(1,101):
        gap=101-2*weight
        prev,value=0,1
        for r,row in enumerate(signed):
            row.append(abs(value))
            if r<T:
                nxt,rem=divmod(gap*value-(100-r)*prev,r+1)
                assert rem==0
                prev,value=value,nxt
    result=prefix(200)[1:T+2]
    for k in range(1,GROUPS+1):
        capacities=[
            3**k*sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
            for w in range(1,101)
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


def anchor_miss(alpha):
    anchor=[comb(ANCHOR_BITS,w)*9**(ANCHOR_BITS-w) if w<=ANCHOR_BITS else 0
            for w in range(T+1)]
    remain=[comb(2*GROUPS*M-ANCHOR_BITS,w)*9**(2*GROUPS*M-ANCHOR_BITS-w) for w in range(T+1)]
    miss=sum(anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
             for d in range(RADIUS+1,min(ANCHOR_BITS,T)+1))
    return miss/DEN


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
        volume += comb(200,t)
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
    if part==1:
        field_checks()
        orbits,fixed=orbit_check()
        w,_,_=derive(10)
        assert w == 663155043702000
        assert w**1000 <= 1<<49237
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
        component=baseline+GROUPS*sum(losses)+anchor+child
        assert 0<=component<=1
        total=1-(1-component)**COMPONENTS
        assert total <= Fraction(953564703246752699,10**18)
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
