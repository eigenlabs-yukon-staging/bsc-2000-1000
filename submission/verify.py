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
    RADIUS,REDUCE_4,BUDGETS,BUDGET_BASE,CAP_BASE,CAP_OFFSETS,PROFILE_BASE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,
    ANCHOR_BITS,PAR_WIDTHS,SYS_WIDTHS,SEARCH_WIDTHS,SYS_ONLY_BASE,
    QR_ROW_BASE,SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from dynamics import derive, data
from joint_admission import child_removed_counts
from intersection import pair_direction_bound
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


def fixed_noise_direction_bound():
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
                cap=[16*3**k*sum((-1)**(k-h)*comb(k,h)*(comb(M*h+ANCHOR_BITS,w)-comb(M*h,w))
                            for h in range(k+1)) for w in range(1,1001)]
            else:
                cap=[3**k*sum((-1)**(k-h)*comb(k,h)*comb(M*h,w)
                            for h in range(k+1)) for w in range(1,1001)]
            n=M*k+ANCHOR_BITS*special
            inside=binomial(n)
            outside=binomial(1000-n)
            con=[0]*(T+1)
            for r,row in enumerate(masks):
                grouped=ranked_groups((v,c) for v,c in zip(row,cap,strict=True) if c)
                for a in range(min(n,T-r)+1):
                    con[a+r]+=envelope(grouped,(7 if special else 1)*inside[a])
            for t in range(T+1):
                result[t]+=comb(GROUPS-1,k)*sum(outside[b]*con[t-b] for b in range(min(t,1000-n)+1))
    assert all(a<=b for a,b in zip(result,result[1:]))
    return result

def pair_rows():
    data=json.loads((HERE/"pairs.json").read_text())
    cutoffs=range(217,227)
    assert set(data)=={str(t) for t in cutoffs}
    return {t:int(data[str(t)]) for t in cutoffs}


@lru_cache(None)
def full_direction_bound():
    rows=pair_rows()
    output=[]
    running=0
    for t,value in enumerate(fixed_noise_direction_bound()):
        if t in rows:
            value=min(value,rows[t])
        running=max(running,value)
        output.append(running)
    return output


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

def anchor_counts():
    complement=binomial(2000-ANCHOR_BITS)
    return [sum(comb(ANCHOR_BITS,d)*complement[t-d]
                for d in range(RADIUS+1,min(t,ANCHOR_BITS)+1))
            for t in range(T+1)]

def root_up(n, d):
    scaled = n << 128
    s = isqrt(scaled//d)
    return s + (s*s*d < scaled)

@lru_cache(None)
def paired_parts(kind,a,b):
    [(w0,qr),(w1,golay),(w2,hamming),(w3,repetition)]=data()
    common=((48,qr,13),(24,golay,1),(8,hamming,2))
    extras=(((1,{(1,0):1,(0,1):1},2),),((1,{(1,1):1},2),),((2,repetition,2),))
    shape=common+extras[kind]
    parts=[]
    for bits,hist,count in shape:
        row=[0]*(bits+1)
        for profile,multiplicity in hist.items():
            h=sum(n*a**(2*i)*b**(2*bits-2*i) for i,n in enumerate(profile))
            root=root_up(h,b**(2*bits))
            for i,n in enumerate(profile):
                row[i]+=multiplicity*n*root
        parts.append((row,count))
    return product(parts,T,4000)


def run_coefficients(length, kind, parity_bits, free, includes_anchor, a, b):
    paired=paired_parts(kind,a,b)
    row=product(((paired,1),(binomial(free),1)),T,4000)
    half_free=(free+1)//2
    scale_num=(a*a+b*b)**half_free
    scale_den=b**(2*half_free)*(1<<((parity_bits-length-2*includes_anchor)//2+18*64))
    return [Fraction(n*scale_num*b**w,scale_den*a**w) for w,n in enumerate(row)]


def geometry(start,length):
    fields=[(start+i)%GROUPS for i in range(length+1)]
    parity=sum(PAR_WIDTHS[i] for i in fields[1:])
    free=parity+sum(SYS_WIDTHS[i] for i in fields[1:-1])
    local=free+SYS_WIDTHS[start]+PAR_WIDTHS[start]
    return min(start,2),parity,free,local,0 in fields[:-1]


def proper_counts(data):
    totals=[Fraction(0)]*(T+1)
    for length in range(1,GROUPS):
        shapes=Counter(geometry(i,length) for i in range(GROUPS))
        for (kind,parity,free,local,includes_anchor),mult in shapes.items():
            by_weight=[Fraction(comb(local,w)) for w in range(T+1)]
            for a in range(28,55,2):
                tilt=run_coefficients(length,kind,parity,free,includes_anchor,a,100)
                by_weight=[min(x,y) for x,y in zip(by_weight,tilt,strict=True)]
            outside=binomial(2000-local)
            for t in range(T+1):
                totals[t]+=mult*sum(x*outside[t-w] for w,x in enumerate(by_weight[:t+1]))
    return [upper_decimal(value,1) for value in totals]


def block_error(child,shapes):
    anchor=anchor_counts()
    proper=proper_counts(shapes)
    last=removed=volume=success=0
    for t,bound in enumerate(full_direction_bound()):
        shell=comb(2000,t)
        volume+=shell
        removed+=min(shell,anchor[t]+child[t]+proper[t])
        available=volume-removed
        if available:
            pairs=volume*bound
            k=1+pairs//(Q*available)
            value=(2*k*Q*available-pairs)//(Q*k*(k+1))
        else:
            value=0
        current=max(last,value)
        assert 0<=current<=Q
        success+=(current-last)*WEIGHTS[t]
        last=current
    return 1-Fraction(success,DEN)


def field_checks():
    x=2
    for j in range(1,M+1):
        x=square(x)
        if j in (2,167):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert sum(SYS_WIDTHS)==sum(PAR_WIDTHS)==1000
    assert len(SEARCH_WIDTHS)==TOTAL_CELLS==36
    assert QR_ROW_BASE+ORBITS*ROW_CAP==BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[12]+1023)==PARITY[12](1023)
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    for e in (0,1):
        assert table_word(SYS_ONLY_BASE+e) == e|(e<<12)|(cost(10,9**(1-e))<<17)
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg and encode(msg).bit_length()<=2000


def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    shapes=data()
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        field_checks()
        orbits,fixed=orbit_check()
        w,_,_=derive(10)
        assert w == 3160110756777792829726456158289815335633316720255983875772330598236424830974228163195405920592926
        assert w**500 <= 1<<160283
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values()),
                          "W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
        return
    if part==1:
        rows=pair_rows()
        for t,declared in rows.items():
            assert pair_direction_bound(t)==declared
        print(json.dumps({"checked_cutoffs":list(rows)}))
        return
    path=json.loads((HERE/"misses.json").read_text())
    assert isinstance(path,list) and len(path)==T+1
    child=[int(value) for value in path]
    assert all(0<=value<=comb(2000,t) for t,value in enumerate(child))
    if part==2:
        assert child_removed_counts(T)==child
        print(json.dumps({"checked_weights":len(child)}))
    else:
        error=block_error(child,shapes)
        assert error<=Fraction(6896914790026339,125000000000000000)
        print(json.dumps({"B_rounded":upper_decimal(error,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
