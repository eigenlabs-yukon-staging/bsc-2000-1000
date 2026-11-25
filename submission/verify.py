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
from construction import (
    COMPONENT_GROUPS, FIRST_WORK_RAM, LIMIT, M, MEMORY_END, MODULUS, PARITY,
    REDUCE_4, SPECS, QR_ROW_END, candidate_row, encode, poly_remainder, square,
    table_word, WIDTHS, WORDS, QR_ROW_BASE, SQUARE_BASE, LIST_BASE,
    QR_PARITY_BASE, SAVED_B, FIELD_COUNT, cost,
)
from dynamics import derive, data, kept_counts
from qr import ORBITS, ROW_CAP, orbit_check, profiles as qr_compute

OUTER_SETUP_BOUND=50000


@lru_cache(None)
def binomial(n,T):
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


def full_direction_bound(spec):
    T=spec.cutoff
    payload=spec.payload
    groups=spec.groups
    signed=[[] for _ in range(T+1)]
    for weight in range(1,payload+1):
        gap=payload+1-2*weight
        prev,value=0,1
        for r,row in enumerate(signed):
            row.append(abs(value))
            if r<T:
                nxt,rem=divmod(gap*value-(payload-r)*prev,r+1)
                assert rem==0
                prev,value=value,nxt
    result=prefix(spec.length)[1:T+2]
    for k in range(1,groups+1):
        capacities=[
            3**k*sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
            for w in range(1,payload+1)
        ]
        assert all(c>=0 for c in capacities)
        sys=binomial(k*M,T)
        contributions=[0]*(T+1)
        for r,row in enumerate(signed):
            table=ranked_groups((v,c) for v,c in zip(row,capacities,strict=True) if c)
            for a in range(min(k*M,T-r)+1):
                contributions[a+r]+=envelope(table,sys[a])
        outside=binomial((groups-k)*M,T)
        for t in range(T+1):
            result[t]+=comb(groups,k)*sum(
                outside[b]*contributions[t-b] for b in range(min(t,(groups-k)*M)+1)
            )
    assert all(a<=b for a,b in zip(result,result[1:]))
    return result


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
def kernel(length,T):
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


def child_shells(spec):
    width_digit=spec.length-2*M+1
    digit_mask=(1<<width_digit)-1
    cutoff=spec.cutoff
    modulus=(1<<(width_digit*(cutoff+1)))-1
    _,caps,_=derive(spec,10,True)
    states=[1]*(LIMIT+1)
    profiles=dict(data())
    for index,width in reversed(list(enumerate(spec.schedule))):
        reserve=len(spec.schedule)-1-index
        next_states=[0]*(LIMIT+1)
        weights=[9**(2*width-i) for i in range(2*width+1)]
        for i,(profile,count) in enumerate(profiles[width].items()):
            S=sum(x*y for x,y in zip(profile,weights,strict=True))
            monomials=Counter()
            for w,amount in enumerate(kept_counts(width,profile)):
                if amount:
                    monomials[cost(S,weights[w]),w]+=amount*count
            sums=[0]*(LIMIT+1)
            for (c,w),multiplicity in monomials.items():
                left=w*width_digit
                for cap in range(c+reserve,LIMIT+1):
                    sums[cap]+=multiplicity*states[cap-c]<<left
            for r,cap in enumerate(caps[index][i]):
                next_states[r]+=sums[cap]
        states=[state&modulus for state in next_states]
    accepted=[[states[cap]>>(w*width_digit)&digit_mask for w in range(cutoff+1)]
              for cap in spec.budgets]
    for row in accepted:
        assert all(x<=comb(spec.length-2*M,w) for w,x in enumerate(row))
    return accepted

def proper_shells(spec):
    T=spec.cutoff
    answer=[Fraction(0)]*(T+1)
    for length in range(1,spec.groups):
        free_bits=M*(2*length-1)
        outside_bits=spec.length-2*M-free_bits
        free=[comb(free_bits,w) for w in range(T+1)]
        outside=[comb(outside_bits,w) for w in range(T+1)]
        directions=kernel(length,T)
        denom=1<<(M*length+1)
        bad=[0]*(T+1)
        for profile,multiplicity in pair_profiles().items():
            for w in range(T+1):
                terms=range(min(w,len(profile)-1)+1)
                all_noise=sum(profile[d]*free[w-d] for d in terms)
                possible=sum(profile[d]*directions[w-d] for d in terms)
                same=sum(profile[d]*free[w-d]*directions[w-d] for d in terms)
                bad[w]+=multiplicity*min(denom*all_noise,all_noise*possible-same)
        for t,value in enumerate(conv(bad,outside)[:T+1]):
            answer[t]+=Fraction(spec.groups*value,denom)
    return answer

def cumulative_error(spec,accepted):
    T=spec.cutoff
    n=spec.length
    Q=1<<spec.payload
    bad=proper_shells(spec)
    anchor=[comb(M,w) for w in range(T+1)]
    child=[comb(n-2*M,w) for w in range(T+1)]
    for d in range(T+1):
        missed=[x-accepted[d][z] if d<=spec.radius else x for z,x in enumerate(child)]
        for t,val in enumerate(conv(anchor,missed)[:T-d+1],d):
            bad[t]+=anchor[d]*val
    volume=removed=previous=suc=0
    for t,D in enumerate(full_direction_bound(spec)):
        vol=comb(n,t)
        volume+=vol
        removed+=min(vol,upper_decimal(bad[t],1))
        A=volume-removed
        if A:
            k=1+volume*D//(Q*A)
            occupied=(2*k*Q*A-volume*D)//(Q*k*(k+1))
            assert occupied<=Q
            if occupied>previous:
                suc+=(occupied-previous)*9**(n-t)
                previous=occupied
    return 1-Fraction(suc,10**n)


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
    assert sum(COMPONENT_GROUPS)==FIELD_COUNT==20
    assert 48==SAVED_B<FIRST_WORK_RAM==68<SQUARE_BASE<QR_ROW_BASE
    assert QR_ROW_BASE+ORBITS*ROW_CAP==QR_ROW_END
    previous=QR_ROW_END
    for spec in SPECS:
        assert spec.budget_base==previous
        assert spec.budget_base+len(spec.budgets)==spec.profile24
        assert spec.profile24+ORBITS==spec.profile1
        assert spec.profile1+2==spec.cap_offsets[0]
        assert len(spec.cap_offsets)==len(spec.schedule)+1
        assert min(spec.budgets)>=len(spec.schedule)
        previous=spec.cap_offsets[-1]
        assert all(table_word(spec.budget_base+d)==c for d,c in enumerate(spec.budgets))
    assert previous==MEMORY_END<=1<<32
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


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    if part==0:
        actual=qr_compute()
        assert actual==data()[0][1]
        field_checks()
        orbits,fixed=orbit_check()
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values()),
                          "memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
        return
    if part==1:
        works={}
        for spec in SPECS:
            w,caps,_=derive(spec,10,True)
            assert len(caps)==len(spec.schedule)
            for layer in caps:
                for profile in layer:
                    assert len(profile)==LIMIT+1
                    assert all(r<=c<=LIMIT for r,c in enumerate(profile))
            works[spec.groups]=str(w)
        print(json.dumps({"one_component_work":works}))
        return
    admitted=json.loads((HERE/"admitted.json").read_text())
    assert set(admitted)=={"3","2"}
    for spec in SPECS:
        rows=admitted[str(spec.groups)]
        assert len(rows)==spec.radius+1
        for row in rows:
            assert len(row)==spec.cutoff+1
            assert all(type(v) is int and 0<=v<=comb(spec.length-2*M,w)
                       for w,v in enumerate(row))
    if part==2:
        for spec in SPECS:
            assert child_shells(spec)==admitted[str(spec.groups)]
        print(json.dumps({"verified_weight_shells":{
            spec.groups:(spec.radius+1)*(spec.cutoff+1) for spec in SPECS}}))
        return
    bound,exponent=[
        Fraction(int(c["score"]["numerator"]),int(c["score"]["denominator"]))
        for c in json.loads((HERE/"submission.json").read_text())["claims"]
    ]
    success=Fraction(1)
    W=OUTER_SETUP_BOUND
    details=[]
    for spec in SPECS:
        component=cumulative_error(spec,admitted[str(spec.groups)])
        assert 0<=component<=1
        count=COMPONENT_GROUPS.count(spec.groups)
        success*=(1-component)**count
        work,_,_=derive(spec,10)
        W+=count*work
        details.append({"groups":spec.groups,"count":count,
                        "component_error_ceil":upper_decimal(component,10**18),
                        "work":str(work)})
    assert 1-success<=bound
    assert W**exponent.denominator<=1<<exponent.numerator
    print(json.dumps({"B_ceil":upper_decimal(1-success,10**18),"W":str(W),
                      "L":str(exponent),"components":details}))


if __name__=='__main__':
    main()
