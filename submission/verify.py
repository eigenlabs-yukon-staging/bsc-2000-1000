"""Exact noise and operation certificate for the split decoder."""
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path
import json
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (
    MODULUS, THETA, WIDTH, SUBMASK, M, EXTENSION_LOW, RADIUS,
    SPLIT_COUNT, ORIENTATION_CHARTS, MEMORY_END, base_product,
    base_power, chart_at, codeword, descriptor, direction, operation_bound,
    parts, recover, recovery, reduction_table, remainder, sqrt, square,
)
from minors import WITNESSES

CUT=240
Q=1<<1000
q=1<<500
G=(q-27216000*2**450)**2
C1=tuple(comb(50,i) for i in range(51))
C=tuple(comb(500,i) for i in range(CUT+1))
C2=tuple(comb(1000,i) for i in range(CUT+1))
D=42
BAD=(2*sum(C[:D]))**2


def prefix(values):
    out=[0]
    for v in values:
        out.append(out[-1]+v)
    return out


@lru_cache(1)
def competitor_bounds():
    balls=prefix(C2)[1:]
    bcoeff=[sum(comb(D,j)*comb(1000-2*D,i-2*j)
                for j in range(min(D,i//2)+1)) for i in range(CUT+1)]
    bias=prefix(bcoeff)[1:]
    double=[]
    for t in range(CUT+1):
        v=0
        for i in range(t+1):
            n,a,b=C2[i],balls[t-i],bias[t-i]
            v+=n*(a+b)+min(n,BAD)*(a-b)
        if v>=Q:
            break
        double.append(v)
    partial=[0]
    for t in range(CUT+1):
        v=0
        for i in range(t//3+1):
            for j in range(i,(t-i)//2+1):
                k=t-i-j
                mult=1 if i==k else 3 if i==j or j==k else 6
                v+=mult*C[i]*C[j]
        total=v+partial[-1]
        if total>=q:
            break
        partial.append(total)
    return double,partial[1:]


def powers(coeffs,max_power,capacity):
    out=[{0:1}]
    for _ in range(max_power):
        nxt={}
        for w,c in out[-1].items():
            for u,d in coeffs:
                if w+u>capacity:break
                nxt[w+u]=nxt.get(w+u,0)+c*d
        out.append(nxt)
    return out


@lru_cache(2)
def group_hist(limit,count):
    hist=[[0]*(count*w//10+1) for w in range(limit+1)]
    for pivot in range(limit//(11-count)+1):
        low=powers(tuple(enumerate(C1[:pivot])),count-1,limit)
        high=powers(tuple((i,C1[i]) for i in range(pivot+1,51)),10-count,limit)
        eqpow=[C1[pivot]**e for e in range(11)]
        for k in range(count):
            for equal in range(count-k,11-k):
                h=10-k-equal
                factor=comb(10,k)*comb(10-k,equal)*eqpow[equal]
                for wlow,v in low[k].items():
                    selected=wlow+(count-k)*pivot
                    part=wlow+equal*pivot
                    for whigh,u in high[h].items():
                        total=part+whigh
                        if total<=limit:
                            hist[total][selected]+=factor*v*u
    for w,row in enumerate(hist):
        assert sum(row)==C[w]
    return hist


def pack(row,bits):
    value=0
    for x in reversed(row):
        value=(value<<bits)|x
    return value


@lru_cache(1)
def parity_tables(limit):
    three=group_hist(limit,3)
    bits=1024
    a=[pack(row,bits) for row in three]
    tables=[]
    for p in range(limit+1):
        v=sum(a[c]*a[p-c] for c in range(p+1))
        mask=(1<<bits)-1
        row=[]
        while v:
            row.append(v&mask)
            v>>=bits
        tables.append(prefix(row))
        assert tables[-1][-1]==C2[p]
    return tables


def certificate(radius):
    double,partial=competitor_bounds()
    stop=len(double)-1
    four=group_hist(stop,4)
    parity=parity_tables(stop)
    tails=[tuple(sum(row[i:]) for i in range(len(row)+1)) for row in four]
    success=[0]*(stop+1)
    def tail(w,l):
        a=tails[w]
        return a[min(max(0,l),len(a)-1)]
    for a in range(min(radius,stop)+1):
        for b in range(a,stop-a+1):
            multiplier=1 if a==b else 2
            ks=[]
            prior=0
            for j in range(min(radius+1,a+4*b//10,b+4*a//10),a-1,-1):
                cumulative=tail(a,j-b)*tail(b,j-a)
                if j<=radius and cumulative>prior:
                    ks.append((radius-j,(cumulative-prior)*multiplier))
                prior=cumulative
            for p in range(stop-a-b+1):
                t=a+b+p
                pa,pb=t-a,t-b
                if pa>=len(partial) or pb>=len(partial):
                    continue
                positive=G-double[t]-q*(partial[pa]+partial[pb])
                if positive<=0:
                    continue
                box=parity[p]
                n=sum(v*box[min(r+1,len(box)-1)] for r,v in ks)
                success[t]+=n*positive
    num=sum(x*9**(2000-t) for t,x in enumerate(success))
    return 1-Fraction(num,Q*10**2000), stop, len(partial)-1

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
    assert max(reduction_table(4)).bit_length()<=64
    return checks,extension


def check_charts():
    rng=random.Random(59267)
    checked=updates=0
    cases=(0,101,9528,ORIENTATION_CHARTS-1,
           ORIENTATION_CHARTS,ORIENTATION_CHARTS+101,SPLIT_COUNT-1)
    for ordinal in cases:
        chart=chart_at(ordinal)
        outer,known,cr,dr=chart
        v=WITNESSES[outer]
        u=rng.randrange(1,1<<M)
        word=codeword(v,u) if outer==0 else codeword(u,v)
        assert square(sqrt(word[0]))==word[0]
        key=(tuple(parts(u)[i] for i in known)
             +tuple(parts(word[2])[i] for i in cr)
             +tuple(parts(word[3])[i] for i in dr))
        assert recover(chart,v,key)==word
        d=descriptor(ordinal)
        assert len(d)==17 and d[0]==outer
        assert tuple(sorted(d[1:11]))==tuple(range(10))
        assert d[1:5]==known and d[11:14]==cr and d[14:17]==dr
        checked+=1
        system=recovery(chart,v)
        for i in (0,49,50,199,200,349,350,499):
            e,dc,dd=direction(chart,v,system,i)
            out=list(word)
            out[1-outer]^=e
            out[2]^=dc
            out[3]^=dd
            change=list(key)
            change[i//WIDTH]^=1<<(i%WIDTH)
            assert recover(chart,v,tuple(change))==tuple(out)
            updates+=1
    return checked,updates


def main():
    json.load(sys.stdin)
    d=json.loads((HERE/"submission.json").read_text())
    b,l=[Fraction(int(x["score"]["numerator"]),int(x["score"]["denominator"]))
         for x in d["claims"]]
    fields=check_field()
    n,o,w=operation_bound()
    assert w**l.denominator<=2**l.numerator and MEMORY_END<=1<<32
    error,global_stop,partial_stop=certificate(RADIUS)
    assert error<=b
    cases,steps=check_charts()
    print(json.dumps({"B":str(b),"L":str(l),"W":str(w),"base_and_extension_gcd":fields,
        "error_ceiling_21_digits":-(-(error.numerator*10**21)//error.denominator),
        "global_end":global_stop,"partial_end":partial_stop,"split_chart_count":SPLIT_COUNT,
        "recovery_cases":cases,"direction_cases":steps,"memory_end_exclusive":MEMORY_END,
        "inner_nodes":str(n),"outer_nodes":str(o)}))


if __name__=="__main__":
    main()
