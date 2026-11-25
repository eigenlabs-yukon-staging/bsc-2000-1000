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
    ANCHOR_BITS,QR_ROW_WORDS,QR_REP_BASE,QR_COST_BASE,poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE)

from dynamics import derive, data
from accelerator import counter_samples
from prefix_cube import cache_samples
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


def pair_rows():
    data=json.loads((HERE/"pairs.json").read_text())
    cutoffs=range(212,231)
    assert set(data)=={str(t) for t in cutoffs}
    return {t:int(data[str(t)]) for t in cutoffs}


@lru_cache(None)
def full_direction_bound():
    legacy=fixed_noise_direction_bound()
    rows=pair_rows()
    result=[]
    current=0
    for t,previous in enumerate(legacy):
        value=min(previous,rows[t]) if t in rows else previous
        current=max(current,value)
        result.append(current)
    return result


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

def proper_counts(data):
    length,a,b=1,42,100
    assert 0 < a < b and 2 * a * a <= b * b
    free = M*(2*length-1)
    local = free+2*M
    parts = []
    for (width, hist), count in zip(data,(20,1,2),strict=True):
        row = [0]*(2*width+1)
        for profile, multiplicity in hist.items():
            h = sum(n*a**(2*i)*b**(4*width-2*i) for i,n in enumerate(profile))
            root = root_up(h,b**(4*width))
            for i,n in enumerate(profile):
                row[i] += multiplicity*n*root
        parts.append((row,count))
    parts.append(([comb(free,i) for i in range(T+1)],1))
    row = product(parts,T,5000)
    weighted = [n*b**i*a**(T-i) for i,n in enumerate(row)]
    outside = [comb(2000-local,i) for i in range(T+1)]
    shells = product(((weighted,1),(outside,1)),T,7000)
    scale = Fraction((a*a+b*b)**(free//2), b**free * (1<<((M*length)//2+23*64)) * a**T)
    return [upper_decimal(GROUPS*scale*n,1) for n in shells]

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
    for i in range(1,M+1):
        x=square(x)
        if i in (100,250):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert QR_ROW_BASE+ORBITS*QR_ROW_WORDS==QR_REP_BASE
    assert QR_ROW_WORDS*2==ROW_CAP==N
    assert QR_REP_BASE+ORBITS==QR_COST_BASE
    assert QR_COST_BASE+14*49==BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((4,0),(4,15),(12,0),(12,2026),(12,4095)):
        row = candidate_row(width,syndrome)
        parts = [(word>>17,(word>>12)&31,word&4095) for word in row]
        assert all(cost > 0 for cost, _, _ in parts)
        assert parts == sorted(parts)
        assert sorted(e for _,_,e in parts) == list(range(1<<width))
        for _,d,e in parts:
            assert d == e.bit_count()+(PARITY[width](e)^syndrome).bit_count()
    msg=(1<<1000)-37
    assert encode(msg)&(Q-1)==msg
    assert encode(msg).bit_length()<=2000


def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    shapes=data()
    rows=pair_rows()
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        field_checks()
        orbits,fixed=orbit_check()
        cube_data = counter_samples()
        cache_data = cache_samples()
        w = int(cache_data["W"])
        assert w == 12090556300019517937004969428759815635665987178972565222314628033298170323706619569510172111815128586782289195893137422017941174208690
        assert w**9677 <= 1 << 4278108
        assert cube_data["cube_ops"] == 6108
        assert (cache_data["setup_ops"], cache_data["reuse_ops"]) == (10273, 2348)
        print(json.dumps({"syndrome_types": len(actual), "syndromes": sum(actual.values()),
                          "orbits": orbits, "fixed": fixed, **cube_data, **cache_data}))
        return
    if part in (1,2):
        cutoffs=range(212,225) if part==1 else range(225,231)
        for t in cutoffs:
            assert pair_direction_bound(t)==rows[t]
        if part==1:
            print(json.dumps({"checked_pair_cutoffs":list(cutoffs)}))
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
        assert error<=Fraction(13793812699701584233989266086719410691211,250000000000000000000000000000000000000000)
        print(json.dumps({"B_rounded":upper_decimal(error,10**42),"scale":10**42}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
