"""Exact eight-chart certificate."""
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path
from math import comb
import random
import sys

HERE=Path('/snapshot/files') if Path('/snapshot/files/construction.py').is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (
    FALLBACK, M, HALF, MODULUS, XI, THETA, LAMBDA, THIRD, RADIUS, CHARTS,
    MEMORY_END, W_bound, base_product, codeword, cube, derivative_columns,
    derivative_entry, directions, good, parts, power, product,
    remainder, split_recover, square,
)

Q=1<<1000
G=((1<<250)-1)**2*((1<<250)-2)**2
CUT=240
BALANCE=42
C=[comb(500,i) for i in range(CUT+1)]
T=[comb(250,i) for i in range(CUT+1)]
D=[comb(1000,i) for i in range(CUT+1)]
BAD=15*(2*sum(C[:BALANCE]))**2


def prefix(values):
    result=[0]
    for x in values:result.append(result[-1]+x)
    return result


@lru_cache(1)
def profiles():
    part=[]
    for s in range(CUT+1):
        n=[0]*(s//2+1)
        for i in range(s+1):
            n[min(i,s-i)]+=T[i]*T[s-i]
        part.append(n)
    tail=[tuple(sum(n[j:]) for j in range(len(n)+1)) for n in part]
    def t(w,u):
        row=tail[w]
        return row[min(max(0,u),len(row)-1)]
    system=[]
    parity=[]
    for w in range(CUT+1):
        sys=[0]*(w+1)
        for a in range(w//2+1):
            b=w-a
            old=0
            for u in range(min(a+b//2,b+a//2),a-1,-1):
                now=t(a,u-b)*t(b,u-a)
                sys[u]+=(now-old)*(1 if a==b else 2)
                old=now
        assert sum(sys)==D[w]
        system.append(tuple((i,n) for i,n in enumerate(sys) if n))
        par=[0]*(w//2+1)
        for a in range(w+1):
            q=C[w-a]
            for i,n in enumerate(part[a]):
                par[i]+=n*q
        assert sum(par)==D[w]
        parity.append(prefix(par))
    balls=prefix(D)[1:]
    bias=prefix(sum(comb(BALANCE,j)*comb(1000-2*BALANCE,i-2*j)
                    for j in range(min(BALANCE,i//2)+1)) for i in range(CUT+1))[1:]
    pair=[]
    for t in range(CUT+1):
        total=0
        for j in range(t+1):
            n,a,h=D[j],balls[t-j],bias[t-j]
            total+=n*(a+h)+min(n,BAD)*(a-h)
        pair.append(total)
    return system,parity,pair


def certificate(r):
    sys,par,pair=profiles()
    volume=prev=success=0
    for t in range(CUT+1):
        for a in range(t+1):
            p=par[t-a]
            for i,n in sys[a]:
                if i>r:break
                volume+=n*p[min(r-i+1,len(p)-1)]
        A=pair[t]
        k=1+A//G
        lower=max(prev,volume*G*G//(Q*(G+A)),volume*(2*k*G-A)//(Q*k*(k+1)))
        assert 0<=prev<=lower<=Q
        success+=(lower-prev)*9**(2000-t)
        prev=lower
    return 1-Fraction(success,10**2000)

def gcd(a,b):
    while b:
        a,b=b,remainder(a,b)
    return a


def field_checks():
    checks={}
    x=2
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
    zeta=power(LAMBDA,THIRD)
    assert zeta!=1 and cube(zeta)==1
    assert good(FALLBACK)
    return checks,trace,hex(zeta)


def check_directions():
    rng=random.Random(590733)
    count=0
    for outer,k,r in CHARTS:
        a,b=(rng.getrandbits(M) for _ in range(2))
        assert good(a) and good(b)
        word=codeword(a,b)
        u,v=word[1-outer],word[outer]
        p,q=parts(u)[k],parts(word[2])[r]
        assert split_recover(v,k,r,p,q)==(u,word[2])
        change=tuple(directions(v,k,r))
        alpha=LAMBDA if outer==0 else 1
        for i in (0,249,250,499):
            e,dc=change[i]
            p1=p^(1<<i) if i<250 else p
            q1=q^(1<<(i-250)) if i>=250 else q
            assert split_recover(v,k,r,p1,q1)==(u^e,word[2]^dc)
            cols=tuple(derivative_columns(alpha,e))
            diff=product(alpha,cube(e))
            for h,part in enumerate(parts(u)):
                for bit in range(HALF):
                    if part>>bit&1:
                        diff^=cols[h*256+bit]
            assert diff==product(alpha,cube(u)^cube(u^e))
            assert derivative_entry(alpha,e,63,255)==cols[504]^cols[505]
            count+=1
    return count


def main():
    json.load(sys.stdin)
    manifest=json.loads((HERE/'submission.json').read_text())
    b,l=[Fraction(int(x['score']['numerator']),int(x['score']['denominator']))
         for x in manifest['claims']]
    field,trace,zeta=field_checks()
    n,o,w=W_bound()
    error=certificate(RADIUS)
    assert error<=b and w**l.denominator<=2**l.numerator and MEMORY_END<=1<<32
    print(json.dumps({'field_checks':field,'quadratic_trace':trace,'noncube_witness':zeta,
        'B':str(b),'L':str(l),'W':str(w),'inner_nodes':str(n),'outer_nodes':str(o),
        'bound_ceiling_21_decimals':-(-(error.numerator*10**21)//error.denominator),
        'memory_end':MEMORY_END,'direction_identity_samples':check_directions()}))


if __name__=='__main__':
    main()
