"""Exact admission shells, ranked parity imbalances and run removals."""
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import random
import sys

HERE=Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (
    AFFINE_BASE,BANDS,DEPTH_BASE,GROUPS,M,MASK,MEMORY_END,MODULUS,
    RADIUS,REDUCE_4,SLICE,SLICES,affine_entry,cube,encode,
    operation_bound,poly_remainder,square,square_entry,table_word,
)

T=230
Q=1<<1000
RADIX_BITS=1120
RADIX_BYTES=RADIX_BITS//8
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


def packed_group(min_systematic,limit,stride):
    table=bytearray((limit+1)*stride*RADIX_BYTES)
    for d in range(min(M,limit//2)+1):
        choose_d=comb(M,d)
        prefix_d=prefix(d)
        for h in range(min(M-d,(limit-2*d)//13)+1):
            threshold=max(0,min(d+1,min_systematic-h))
            ways=(1<<d)-prefix_d[threshold]
            value=choose_d*comb(M-d,h)*ways
            offset=((2*d+13*h)*stride+h)*RADIX_BYTES
            table[offset:offset+RADIX_BYTES]=value.to_bytes(RADIX_BYTES,"little")
    return int.from_bytes(table,"little")


def band_shells(index):
    low,high,budget=BANDS[index]
    limit=budget-9*(GROUPS-1)
    assert limit>=0 and limit//2<=750
    assert 4*comb(1500,limit//2)<1<<RADIX_BITS
    stride=1+limit//13
    slots=(limit+1)*stride
    mask=(1<<(RADIX_BITS*slots))-1
    p=packed_group(low,limit,stride)
    q=packed_group(high+1,limit,stride)
    result=(((p+q)*((p*p+q*q)&mask))&mask).to_bytes(slots*RADIX_BYTES,"little")
    inner=[0]*(T+1)
    for d in range(min(limit//2,T)+1):
        for h in range(min((limit-2*d)//13,(T-d)//2)+1):
            offset=((2*d+13*h)*stride+h)*RADIX_BYTES
            inner[d+2*h]+=int.from_bytes(result[offset:offset+RADIX_BYTES],"little")
    outer=[
        sum(comb(M,i)*comb(M,t-i) for i in range(low,min(high,t)+1))
        for t in range(T+1)
    ]
    return [sum(outer[i]*inner[t-i] for i in range(t+1)) for t in range(T+1)]


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


def run_conflict_weights(length):
    local=(2*length+1)*M
    outside=binomial(2000-local)
    free=binomial(local-2*M)
    part=binomial(M)
    coeff=[sum(part[h]*free[w-2*h] for h in range(w//2+1)) for w in range(T+1)]
    masses=[0]*(T+1)
    denominator=1<<(length*M)
    for j in range(min(M,T)+1):
        volume=0
        for w in range(T+1-j):
            count=coeff[w]<<j
            volume+=count
            masses[j+w]+=part[j]*count*min(denominator,volume<<length)
        for i in range(2,T+1):
            coeff[i]-=coeff[i-2]
            assert coeff[i]>=0
    weights=[sum(outside[v-w]*masses[w] for w in range(v+1)) for v in range(T+1)]
    return weights,denominator


def certificate(shells):
    pairs=full_direction_bound()
    assert all(0<=n<=comb(2000,i) for i,n in enumerate(shells))
    previous=volume=success=0
    lower=[]
    alpha=[]
    for count,bound,weight in zip(shells,pairs,WEIGHTS,strict=True):
        volume+=count
        k=1+bound//Q
        current=max(previous,volume*(2*k*Q-bound)//(Q*k*(k+1)))
        assert current<=Q
        success+=(current-previous)*weight
        previous=current
        lower.append(current)
        alpha.append(Fraction(2,k+1))
    loss=Fraction(0)
    details=[]
    for length in range(1,GROUPS):
        weights,q=run_conflict_weights(length)
        term=sum(count*weight*scale for count,weight,scale in zip(weights,WEIGHTS,alpha,strict=True))/(q*DEN)
        details.append(term)
        loss+=GROUPS*term
    return 1-Fraction(success,DEN)+loss,success,details,lower


def field_and_table_checks():
    x=2
    check={}
    for i in range(1,M+1):
        x=square(x)
        if i in (50,125):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
            check[i]=a
    assert x==2
    rng=random.Random(577013)
    for i in (0,17,18,63,64,128,249):
        x=rng.getrandbits(M)
        sq=0
        for k in range((M+7)//8):
            sq^=square_entry(k,(x>>(8*k))&255)
        assert sq==square(x)
        delta=0
        for k in range(SLICES):
            delta^=affine_entry(i,k,(x>>(SLICE*k))&((1<<SLICE)-1))
        assert cube(x)^delta==cube(x^(1<<i))
    assert SLICES*SLICE<=256
    assert 48<REDUCE_4<DEPTH_BASE<AFFINE_BASE<MEMORY_END<=1<<32
    assert max(table_word(REDUCE_4+i) for i in range(16)).bit_length()<=64
    assert 81**2<=1<<13 and 82**(2*M)<=(1<<9)*81**(2*M)
    assert BANDS[0][0]==0 and all(a[1]+1==b[0] for a,b in zip(BANDS,BANDS[1:]))
    assert BANDS[-1][1]==RADIUS<=M
    for lo,hi,budget in BANDS:
        assert budget>=9*(GROUPS-1)
        for depth in range(lo,hi+1):
            assert table_word(DEPTH_BASE+depth)==budget
    msg=(1<<1000)-1
    encoded=encode(msg)
    assert encoded.bit_length()<=2000 and encoded&msg==msg
    return check


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


def main():
    band=json.load(sys.stdin)["band"]
    stored=json.loads((HERE/"shells.json").read_text())
    shells=band_shells(band)
    assert shells==stored[band]
    report={"band":band,"upper_weight":T,"checked_coefficients":len(shells)}
    if band==3:
        field=field_and_table_checks()
        combined=[sum(row) for row in zip(*stored,strict=True)]
        error,numerator,loss,lower=certificate(combined)
        manifest=json.loads((HERE/"submission.json").read_text())
        B,L=[Fraction(int(c["score"]["numerator"]),int(c["score"]["denominator"]))
             for c in manifest["claims"]]
        W,counts=operation_bound()
        assert error<=B and W**L.denominator<=2**L.numerator
        report.update({
            "B":str(B),"L":str(L),"W":str(W),"field":field,"memory_end":MEMORY_END,
            "bands":BANDS,"count_terms":counts,"success_numerator_before_removal":str(numerator),
            "last_occupied_count":str(lower[-1]),
            "one_start_removal_ceil_1e24":[str(upper_decimal(x,10**24)) for x in loss],
        })
    print(json.dumps(report))


if __name__=="__main__":
    main()
