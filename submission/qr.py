from collections import Counter
from array import array
from functools import lru_cache
from math import comb
import sys

GEN = 9205487
D = DIM = 24
N = 1<<D
I = (1<<47)-1
MASK = (1<<31)-1
ORBITS = 426
PERMUTATIONS = 51888
ROW_CAP = 7<<18
MOD = 47
INF = 47
INV = tuple(0 if j==0 else pow(j,-1,MOD) for j in range(MOD))
PIVOTS = tuple(range(23))+(47,)

def remainder(x):
    while x.bit_length()>23:
        x^=GEN<<(x.bit_length()-24)
    return x

def parity(x):
    y=remainder(x<<23)
    return y|(((x.bit_count()+y.bit_count())&1)<<23)

COLS=tuple(1<<i for i in range(23))+tuple(parity(1<<i) for i in range(24))+(1<<23,)

def syndrome(word):
    result=0
    while word:
        b=word&-word
        result^=COLS[b.bit_length()-1]
        word^=b
    return result

def vector(a,b,c,d,x):
    if x==INF:
        return a*INV[c]%MOD if c else INF
    den=(c*x+d)%MOD
    return (a*x+b)*INV[den]%MOD if den else INF

def rank(rows):
    pivots=[0]*DIM
    result=0
    for row in rows:
        while row:
            pivot=row.bit_length()-1
            old=pivots[pivot]
            if old:
                row^=old
            else:
                pivots[pivot]=row
                result+=1
                break
    return result

def matrices():
    for a in range(MOD):
        for b in range(MOD):
            if a:
                for c in range(MOD):
                    d=((1+b*c)*INV[a])%MOD
                    mat=(a,b,c,d)
                    if mat<=tuple(-n%MOD for n in mat):
                        yield mat
            elif b:
                c=-INV[b]%MOD
                for d in range(MOD):
                    mat=(a,b,c,d)
                    if mat<=tuple(-n%MOD for n in mat):
                        yield mat

def orbit_check():
    code=[(GEN<<i)|(1<<47) for i in range(DIM)]
    for g in ((1,1,0,1),(0,46,1,0)):
        p=[vector(*g,i) for i in range(48)]
        assert sorted(p)==list(range(48))
        assert all(syndrome(sum(1<<p[i] for i in range(48) if word>>i&1))==0 for word in code)
    totals=Counter()
    for mat in matrices():
        dim=DIM-rank((COLS[vector(*mat,x)]^(1<<j) for j,x in enumerate(PIVOTS)))
        totals[dim]+=1
    count=sum(totals.values())
    assert count==51888
    n,rem=divmod(sum(freq<<dim for dim,freq in totals.items()),count)
    assert rem==0
    assert n==ORBITS
    return n,dict(sorted(totals.items()))

def quotient(a,b):
    q=0
    while a.bit_length()>=b.bit_length():
        shift=a.bit_length()-b.bit_length()
        a ^= b<<shift
        q |= 1<<shift
    assert a==0
    return q

def linear_values(rows):
    x=[0]
    for v in rows:
        x.extend(v^old for old in x[:])
    return x

def weights():
    rows=[(GEN<<i)|(1<<47) for i in range(D)]
    assert all((a&b).bit_count()%2==0 for a in rows for b in rows)
    lo=linear_values(rows[:12])
    hi=linear_values(rows[12:])
    answer=bytes((a^b).bit_count() for a in hi for b in lo)
    assert set(answer)=={0,12,16,20,24,28,32,36,48}
    return answer

def wht(data,weight,unit):
    records=bytearray(4*N)
    records[::4]=data.translate(bytes(int(i==weight) for i in range(256)))
    packed=int.from_bytes(records,'little')
    del records
    for j in range(D):
        half=1<<j
        shift=32*half
        mask=int.from_bytes((b'\xff\xff\xff\x7f'*half+b'\0'*(4*half))*(N//(2*half)),'little')
        lo=packed&mask
        hi=(packed>>shift)&mask
        del packed
        guard=(unit&mask)<<31
        plus=(lo+hi)&mask
        minus=(lo+guard-hi)&mask
        del lo,hi,guard,mask
        packed=plus|(minus<<shift)
        del plus,minus
    result=array('I')
    assert result.itemsize==4
    result.frombytes(packed.to_bytes(4*N,'little'))
    if sys.byteorder!='little':
        result.byteswap()
    return result

def profiles():
    data=weights()
    parity_label=quotient(I,GEN)
    assert parity_label.bit_length()<=D
    unit=int.from_bytes(b'\x01\x00\x00\x00'*N,'little')
    transforms=[]
    for w in (12,16,20,24):
        transforms.append(wht(data,w,unit))
    del data,unit
    counts=Counter((a,b,c,d,(i&parity_label).bit_count()&1)
        for i,(a,b,c,d) in enumerate(zip(*transforms,strict=True)))
    values=(0,12,16,20,24,28,32,36,48)
    basis=[[sum((-1)**j*comb(w,j)*comb(48-w,d-j) for j in range(min(w,d)+1)) for d in range(49)]
           for w in values]
    profiles=Counter()
    for (a,b,c,d,odd),mult in counts.items():
        a,b,c,d=(v if v<1<<30 else v-(1<<31) for v in (a,b,c,d))
        sig=-1 if odd else 1
        row=[]
        for h in range(49):
            value=sum(v*p[h] for v,p in zip((1,a,b,c,d,sig*c,sig*b,sig*a,sig),basis,strict=True))
            q,r=divmod(value,N)
            assert q>=0 and r==0
            row.append(q)
        profiles[tuple(row)] += mult
    assert len(profiles)==14
    assert sum(profiles.values())==N
    return profiles

@lru_cache(None)
def group():
    return tuple(matrices())


@lru_cache(None)
def permutations():
    return tuple(tuple(vector(*matrix,x) for x in range(48)) for matrix in group())


@lru_cache(None)
def orbit_labels():
    labels=array('I',[0xffffffff])*N
    representatives=[]
    perm=permutations()
    columns=tuple(tuple(COLS[p[j]] for j in PIVOTS) for p in perm)
    for r in range(N):
        if labels[r] != 0xffffffff:
            continue
        index=len(representatives)
        representatives.append(r)
        active=[i for i in range(D) if r>>i&1]
        for k,cols in enumerate(columns):
            s=0
            for i in active:
                s ^= cols[i]
            if labels[s]==0xffffffff:
                labels[s]=index|(k<<9)
    assert len(representatives)==ORBITS
    return labels,tuple(representatives)


def permute_byte(k,address,value):
    p=permutations()[k]
    return sum(1<<p[8*address+i] for i in range(8) if value>>i&1)


def error_pattern(e,s):
    p=parity(e)^s
    return (p&((1<<23)-1))|(e<<23)|((p>>23)<<47)


@lru_cache(maxsize=1)
def row_words(index,precision):
    representative=orbit_labels()[1][index]
    buckets=[array('Q') for _ in range(49)]
    for e in range(N):
        val=error_pattern(e,representative)
        buckets[val.bit_count()].append(val)
    for bucket in buckets:
        ordered=sorted(bucket)
        bucket[:]=array('Q',ordered)
    s=sum(len(row)*9**(48-d) for d,row in enumerate(buckets))
    output=array('Q')
    for d,row in enumerate(buckets):
        cost=((s**precision-1)//9**((48-d)*precision)).bit_length()
        take=min(ROW_CAP-len(output),len(row))
        output.extend(e|(cost<<48) for e in row[:take])
        if len(output)==ROW_CAP:
            break
    assert len(output)==ROW_CAP
    return output
