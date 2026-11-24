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
    AS_map, AS_pivots, MODULUS, LAMBDA, THETA, M, WIDTH, SUBMASK, EXTENSION_LOW, BYTES_PER_CHUNK,
    MEMORY_END, CHART_COSTS, PAIRS, RADII, SPLIT_RADIUS, SPLIT_CHARTS,
    THIRD, CUBE_EXP, W_bound, base_product, base_power, chart, codeword,
    constant, cube, derivative_columns, affine_entry, directions, evaluate_product, evaluation_columns, interpolation_columns, linear_slice_entry, nonzero_minor_count,
    inverse, map_tags, packed_words, parts, join, power, product, reduction_table,
    remainder, split_map, split_recover, square,
)

G=(2**500-317520*2**450)**2
Q=2**1000
CUT=240
D=42
C=[comb(500,j) for j in range(CUT+1)]
C1000=[comb(1000,j) for j in range(CUT+1)]
C2000=[comb(2000,j) for j in range(CUT+1)]
H=[comb(50,j) for j in range(51)]


def prefix(seq):
    p=[0]
    for x in seq:p.append(p[-1]+x)
    return p


def powers(coefficients,capacity):
    polys=[{0:1}]
    for _ in range(5):
        nxt={}
        for w,c in polys[-1].items():
            for t,v in coefficients:
                if w+t>capacity:
                    break
                nxt[w+t]=nxt.get(w+t,0)+c*v
        polys.append(nxt)
    return polys


def tables(limit):
    histogram=[[0]*(t//2+1) for t in range(limit+1)]
    for pivot in range(limit//6+1):
        low=powers(tuple(enumerate(H[:pivot])),limit)
        high=powers(tuple((i,H[i]) for i in range(pivot+1,51)),limit)
        eqpow=[H[pivot]**i for i in range(11)]
        for k in range(5):
            for e in range(5-k,11-k):
                h=10-k-e
                factor=comb(10,k)*comb(10-k,e)*eqpow[e]
                baseline=e*pivot
                for a,c in low[k].items():
                    smallest=a+(5-k)*pivot
                    for b,d in high[h].items():
                        total=a+b+baseline
                        if total<=limit:
                            histogram[total][smallest]+=factor*c*d
    tails=[]
    for weight,h in enumerate(histogram):
        suffix=tuple(sum(h[j:]) for j in range(len(h)))+(0,)
        assert suffix[0]==C[weight]
        tails.append(suffix)
    return tails


def certificate(rsys,rmul,rcube,rpar,rsplit):
    tails=tables(CUT-min(rsys,rpar)-1)
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


def trim(poly):
    while poly and not poly[-1]:
        poly.pop()
    return poly


def polynomial_remainder(a,b):
    out=trim(list(a))
    inv=base_power(b[-1],SUBMASK-1)
    while len(out)>=len(b):
        offset=len(out)-len(b)
        factor=base_product(out[-1],inv)
        for i,v in enumerate(b):
            out[offset+i]^=base_product(factor,v)
        trim(out)
    return out


def base_polynomial_gcd(a,b):
    a,b=trim(list(a)),trim(list(b))
    while b:
        a,b=b,polynomial_remainder(a,b)
    leading=base_power(a[-1],SUBMASK-1)
    return [base_product(leading,x) for x in a]


def check_field():
    x=2
    checks={}
    for i in range(1,WIDTH+1):
        x=base_product(x,x)
        if i in (10,25):
            checks[i]=gcd(x^2,MODULUS)
            assert checks[i]==1
    assert x==2
    t=THETA
    coeff=dict(EXTENSION_LOW)
    modulus=[coeff.get(i,0) for i in range(10)]+[1]
    extension={}
    for e in range(1,11):
        for _ in range(WIDTH):
            t=square(t)
        if e in (2,5):
            extension[e]=base_polynomial_gcd(modulus,parts(t^THETA))
            assert extension[e]==[1]
    assert t==THETA
    zeta=constant("zeta")
    assert zeta!=1 and cube(zeta)==1 and zeta==power(LAMBDA,THIRD)
    assert 3*CUBE_EXP==2*THIRD+1
    assert CUBE_EXP==1+56*sum(1<<(6*i) for i in range(83))
    assert product(LAMBDA,constant("lambda_inv"))==1
    assert max(reduction_table(4)).bit_length()<=64
    minors=nonzero_minor_count()
    assert minors==63504
    pivots=AS_pivots()
    assert len(pivots)==M-1
    for _,image,_ in pivots:
        preimage=AS_map(image)
        assert square(preimage)^preimage==image
    return checks,extension,minors


def check_charts():
    rng=random.Random(31749)
    tested=split=updates=0
    for trial in range(3):
        a,b=(rng.randrange(1,1<<M) for _ in range(2))
        word=codeword(a,b)
        for i,j in PAIRS:
            found=chart((i,j),word[i],word[j])
            assert (a,b) in found and len(found)<=3
            assert all(a0 and b0 and (codeword(a0,b0)[i],codeword(a0,b0)[j])==(word[i],word[j]) for a0,b0 in found)
            tested+=1
        for number in (0,263,10019,47335,63503,63504,81807,107100,127007):
            outer,known,rows=SPLIT_CHARTS[number]
            u=word[1-outer]
            v=word[outer]
            key=tuple(parts(u)[i] for i in known)+tuple(parts(word[2])[i] for i in rows)
            assert split_recover(v,known,rows,key)==(u,word[2])
            split+=1
            if trial!=0 or number not in (0,63503,63504,127007):
                continue
            ds=tuple(directions(v,known,rows))
            for i in (0,49,50,249,250,449,450,499):
                delta,c_delta=ds[i]
                key1=list(key)
                key1[i//WIDTH]^=1<<(i%WIDTH)
                assert split_recover(v,known,rows,key1)==(u^delta,word[2]^c_delta)
                alpha=LAMBDA if outer==0 else 1
                cols=tuple(derivative_columns(alpha,delta))
                diff=product(alpha,cube(delta))
                for j,part in enumerate(parts(u)):
                    for bit in range(WIDTH):
                        if part>>bit&1:
                            diff^=cols[j*WIDTH+bit]
                assert diff==product(alpha,cube(u)^cube(u^delta))
                if i==499:
                    expected=0
                    for v0 in cols[489:500]:
                        expected^=v0
                    assert affine_entry(alpha,delta,39,2047)==expected
                assert all(word>>50==0 for word in packed_words(diff))
                updates+=1
    return tested,split,updates


def evaluate_examples():
    rng=random.Random(57031)
    pairs=[(0,1),(1<<499,1<<460)]
    pairs.extend((rng.getrandbits(M),rng.getrandbits(M)) for _ in range(14))
    for a,b in pairs:
        assert evaluate_product(a,b)==product(a,b)
        for columns in (evaluation_columns(),interpolation_columns()):
            lane=rng.randrange(len(columns))
            v=rng.getrandbits(WIDTH)
            chunks=[linear_slice_entry(columns,offset,bits,lane,(v>>offset)&((1<<bits)-1))
                    for offset,bits in ((0,20),(20,20),(40,10))]
            assert tuple(x^y^z for x,y,z in zip(*chunks))==tuple(base_product(v,c) for c in columns[lane])
    return len(pairs)


def main():
    json.load(sys.stdin)
    manifest=json.loads((HERE/'submission.json').read_text())
    B,L=[Fraction(int(c["score"]["numerator"]),int(c["score"]["denominator"]))
         for c in manifest["claims"]]
    checked,extension,minors=check_field()
    n,o,ns,os,w=W_bound()
    b,rejects=certificate(RADII[0],RADII[1],RADII[2],RADII[5],SPLIT_RADIUS)
    assert b<=B and w**L.denominator<=2**L.numerator
    assert MEMORY_END<=1<<32
    tested,split,updates=check_charts()
    print(json.dumps({"base_gcd_checks":checked,"extension_gcd_checks":extension,"minor_witness_count":minors,
        "zeta_hex":hex(constant("zeta")),"error_bound":str(B),
        "error_ceiling_at_21_places":-(-(b.numerator*10**21)//b.denominator),
        "rejected_shells_nonempty":sum(x>0 for x in rejects),
        "W":str(w),"L":str(L),"memory_end":MEMORY_END,
        "field_map_count":len(map_tags()),"ordinary_chart_cases":tested,
        "split_chart_cases":split,"direction_update_cases":updates,"product_examples":evaluate_examples(),
        "ordinary_inner_counts":[str(x) for x in n],"ordinary_outer_counts":[str(x) for x in o],
        "split_inner_count":str(ns),"split_outer_count":str(os),"split_chart_count":len(SPLIT_CHARTS)}))


if __name__=='__main__':
    main()
