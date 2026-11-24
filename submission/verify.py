"""Exact finite certificate for ordinary and split information charts."""
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (
    AS_map, AS_pivots, MODULUS, LAMBDA, XI, THETA, M, HALF, SUBMASK,
    MEMORY_END, CHART_COSTS, PAIRS, RADII, SPLIT_RADIUS, SPLIT_CHARTS,
    THIRD, CUBE_EXP, W_bound, base_product, base_power, chart, codeword,
    constant, cube, derivative_columns, derivative_entry, directions, good,
    inverse, map_tags, packed_words, parts, power, product, reduction_table,
    remainder, split_recover, square,
)

G=((2**250-1)*(2**250-2))**2
Q=2**1000
CUT=240
D=42
C=[comb(500,j) for j in range(CUT+1)]
C1000=[comb(1000,j) for j in range(CUT+1)]
C2000=[comb(2000,j) for j in range(CUT+1)]
H=[comb(250,j) for j in range(CUT+1)]


def prefix(seq):
    p=[0]
    for x in seq:p.append(p[-1]+x)
    return p


def tables():
    tails=[]
    for w in range(CUT+1):
        p=prefix(H[j]*H[w-j] for j in range(w+1))
        tails.append(tuple(p[w-k+1]-p[k] for k in range(w//2+1))+(0,))
    return tails


def certificate(rsys,rmul,rcube,rpar,rsplit):
    tails=tables()
    amax=CUT-rpar-1
    pmax=CUT-rsys-1
    rejects=[0]*(CUT+1)

    def tail(w,k):
        series=tails[w]
        return series[min(max(0,k),len(series)-1)]

    @lru_cache(None)
    def parity(a,k):
        low=max(0,rmul-a+1)
        high_d=max(0,rcube-a+1)
        return tuple(sum(tail(c,k)*C[p-c] for c in range(low,p-high_d+1))
                     for p in range(rpar+1,pmax+1))

    for s in range(rsys+1,amax+1):
        max_idx=CUT-s-rpar
        for a in range(s//2+1):
            b=s-a
            previous=0
            scale=1 if a==b else 2
            for k in range(min(a+b//2,b+a//2),a-1,-1):
                surviving=tail(a,k-b)*tail(b,k-a)
                n=(surviving-previous)*scale
                previous=surviving
                if n:
                    for i,v in enumerate(parity(a,max(0,rsplit-k+1))[:max_idx]):
                        rejects[s+i+rpar+1]+=n*v
            assert previous==C[a]*C[b]

    ordinary=prefix(C1000)[1:]
    coefficient=[sum(comb(D,j)*comb(1000-2*D,k-2*j)
                     for j in range(min(D,k//2)+1)) for k in range(CUT+1)]
    bias=prefix(coefficient)[1:]
    bad=15*(2*sum(C[:D]))**2
    V=0
    prior=0
    success=0
    for t in range(CUT+1):
        V+=C2000[t]-rejects[t]
        A=0
        for s in range(t+1):
            n=C1000[s]
            v,b=ordinary[t-s],bias[t-s]
            A+=n*(v+b)+min(n,bad)*(v-b)
        k=1+A//G
        low=max(prior,V*G*G//(Q*(G+A)),V*(2*k*G-A)//(Q*k*(k+1)))
        assert prior <= low <= Q
        success+=(low-prior)*9**(2000-t)
        prior=low
    return 1-Fraction(success,10**2000), rejects

def gcd(a, b):
    while b:
        a,b=b,remainder(a,b)
    return a


def check_field():
    x=2
    checks={}
    for i in range(1,HALF+1):
        x=base_product(x,x)
        if i in (50,125):
            checks[i]=gcd(x^2,MODULUS)
            assert checks[i]==1
    assert x==2
    x=XI
    trace=0
    for _ in range(HALF):
        trace^=x
        x=base_product(x,x)
    assert trace==1 and square(THETA)^THETA==XI
    zeta=constant("zeta")
    assert zeta!=1 and cube(zeta)==1 and zeta==power(LAMBDA,THIRD)
    assert 3*CUBE_EXP==2*THIRD+1
    assert CUBE_EXP==1+56*sum(1<<(6*i) for i in range(83))
    assert product(LAMBDA,constant("lambda_inv"))==1
    assert max(reduction_table(4)).bit_length()<=64
    pivots=AS_pivots()
    assert len(pivots)==M-1
    for _,image,_ in pivots:
        x=AS_map(image)
        assert square(x)^x==image
    return checks,trace


def check_charts():
    rng=random.Random(266310)
    tested=0
    split=0
    updates=0
    for _ in range(4):
        a,b=(rng.getrandbits(M) for _ in range(2))
        assert good(a) and good(b)
        word=codeword(a,b)
        for i,j in PAIRS:
            found=chart((i,j),word[i],word[j])
            assert (a,b) in found and len(found)<=3
            assert all(good(a0) and good(b0) and (codeword(a0,b0)[i],codeword(a0,b0)[j])==(word[i],word[j]) for a0,b0 in found)
            tested+=1
        for outer,j,k in SPLIT_CHARTS:
            u=word[1-outer]
            v=word[outer]
            p,q=parts(u)[j],parts(word[2])[k]
            assert split_recover(v,j,k,p,q)==(u,word[2])
            split+=1
            if _ != 0:
                continue
            ds=tuple(directions(v,j,k))
            for i in (0,249,250,499):
                delta,c_delta=ds[i]
                p1=p^(1<<i) if i<250 else p
                q1=q^(1<<(i-250)) if i>=250 else q
                assert split_recover(v,j,k,p1,q1)==(u^delta,word[2]^c_delta)
                alpha=LAMBDA if outer==0 else 1
                cols=tuple(derivative_columns(alpha,delta))
                diff=product(alpha,cube(delta))
                for h,part in enumerate(parts(u)):
                    for bit in range(250):
                        if part>>bit&1:
                            diff^=cols[h*256+bit]
                assert diff==product(alpha,cube(u)^cube(u^delta))
                assert derivative_entry(alpha,delta,63,255)==cols[504]^cols[505]
                packed=tuple(packed_words(diff))
                assert packed[3]>>58==packed[7]>>58==0
                updates+=1
    return tested,split,updates


def main():
    json.load(sys.stdin)
    manifest=json.loads((HERE/'submission.json').read_text())
    B,L=[Fraction(int(c["score"]["numerator"]),int(c["score"]["denominator"]))
         for c in manifest["claims"]]
    checked,trace=check_field()
    n,o,w=W_bound()
    b,rejects=certificate(RADII[0],RADII[1],RADII[2],RADII[5],SPLIT_RADIUS)
    assert b<=B and w**L.denominator<=2**L.numerator
    assert MEMORY_END<=1<<32
    tested,split,updates=check_charts()
    print(json.dumps({"base_gcd_checks":checked,"quadratic_trace":trace,
        "zeta_hex":hex(constant("zeta")),"error_bound":str(B),
        "error_ceiling_at_21_places":-(-(b.numerator*10**21)//b.denominator),
        "rejected_shells_nonempty":sum(x>0 for x in rejects),
        "W":str(w),"L":str(L),"memory_end":MEMORY_END,
        "field_map_count":len(map_tags()),"ordinary_chart_cases":tested,
        "split_chart_cases":split,"direction_update_cases":updates,
        "inner_node_counts":[str(x) for x in n],"outer_node_counts":[str(x) for x in o]}))


if __name__=='__main__':
    main()
