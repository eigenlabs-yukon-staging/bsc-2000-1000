"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import json
from math import comb, factorial, prod
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import (GROUPS,LIMIT,M,MEMORY_END,MODULUS,PARITY,PRECISION,
    RADIUS,REDUCE_4,SLOPE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS)

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


def envelope_sum(table,capacity):
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
                contributions[a+r]+=envelope_sum(table,sys[a])
        outside=binomial((GROUPS-k)*M)
        for t in range(T+1):
            result[t]+=comb(GROUPS,k)*sum(
                outside[b]*contributions[t-b] for b in range(min(t,(GROUPS-k)*M)+1)
            )
    assert all(a<=b for a,b in zip(result,result[1:]))
    return result


def profiles(width, parity_map):
    weights = [i.bit_count() for i in range(1 << width)]
    values = [parity_map(i) for i in range(1 << width)]
    histogram = Counter()
    for syndrome in range(1 << width):
        profile = [0] * (2*width+1)
        for w, p in zip(weights, values):
            profile[w + weights[p ^ syndrome]] += 1
        histogram[tuple(profile)] += 1
    return histogram

DATA = [(w, profiles(w,PARITY[w])) for w in (12,4)]

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

def admission(precision, limit, radius):
    parts = []
    for (width, hist), exponent in zip(DATA,(72,9)):
        row = [0]*(limit+1)
        weights = [9**(2*width-d) for d in range(2*width+1)]
        for profile, multiplicity in hist.items():
            normalizer = sum(a*b for a,b in zip(profile,weights))
            for n, weight in zip(profile,weights):
                if n:
                    c = cost(normalizer,weight)
                    assert precision == PRECISION and c > 0
                    if c <= limit:
                        row[c] += multiplicity*n*weight
        parts.append((row,exponent))
    row = product(parts,limit,6000)
    for i in range(1,len(row)):
        row[i] += row[i-1]
    admitted = sum(comb(M,i)*9**(M-i)*row[limit-SLOPE*i] for i in range(radius+1))
    return 1-Fraction(admitted,10**1900)

def splits(n, k):
    if k == 1:
        yield (n,)
    else:
        for s in range(n+1):
            for tail in splits(n-s,k-1):
                yield (s,)+tail

def packed(row, stride):
    return int.from_bytes(b''.join(x.to_bytes(stride,'little') for x in row),'little')


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
def direction_bound(m,a):
    p=m*a
    ordinary=(a-1)*m
    result=[sum(binomial(p+ordinary)[:t+1]) for t in range(T+1)]
    rows=[[] for _ in range(T+1)]
    for w in range(1,p+1):
        g=p+1-2*w
        prev,cur=0,1
        for l,row in enumerate(rows):
            row.append(abs(cur))
            if l<T:
                nxt,rem=divmod(g*cur-(p-l)*prev,l+1)
                assert rem==0
                prev,cur=cur,nxt
    for k in range(a):
        outside=binomial(ordinary-k*m)
        degree=k*m
        multi=[(3**k)*sum((-1)**(k-h)*comb(k,h)*comb(h*m,w)
                          for h in range(1,k+1)) for w in range(1,degree+1)]
        syst=binomial(degree)
        total=[0]*(T+1)
        for l,row in enumerate(rows):
            if k==0:
                total[l]+=max(row[:m])
                continue
            shifted=shifted_max(tuple(row),m)
            # shifted[S] runs over exactly S+1 through S+m.
            entries_nonfirst=ranked_groups((row[w-1],multi[w-1]) for w in range(1,degree+1) if multi[w-1])
            entries_withfirst=ranked_groups((shifted[w],multi[w-1]) for w in range(1,degree+1) if multi[w-1])
            for x in range(min(T-l,degree)+1):
                cap=syst[x]
                total[l+x]+=envelope_sum(entries_nonfirst,cap)+envelope_sum(entries_withfirst,cap)
        choose=comb(a-1,k)
        for n in range(T+1):
            result[n]+=choose*sum(outside[b]*total[n-b] for b in range(min(n,ordinary-degree)+1))
    return result


def loss(a):
    m=M
    local=(2*a+1)*m
    free=(2*a-1)*m
    threshold=1<<(m*a+1)
    outside=[comb(2000-local,i) for i in range(T+1)]
    unpaired=[comb(free,i) for i in range(T+1)]
    correction=direction_bound(m,a)
    envelope = [n+(correction[i-1] if i else 0) for i,n in enumerate(correction)]
    capacity=(1<<m)*max(sum(unpaired),sum(envelope))
    stride=(capacity.bit_length()+7)//8
    nbytes=(T//2+1)*stride
    modmask=(1<<(8*nbytes))-1
    free_even=packed(unpaired[::2],stride)
    free_odd=packed(unpaired[1::2],stride)
    rank_even=packed(envelope[::2],stride)
    rank_odd=packed(envelope[1::2],stride)
    totals=[0]*(T+1)
    groups=[]
    for (_, hist), count in zip(DATA,(8,1)):
        types=[]
        for row,mult in hist.items():
            d=next(i for i,n in enumerate(row) if n)
            assert all(not row[i] for i in range(d+1,len(row),2))
            types.append((d,packed(row[d::2],stride),mult))
        powers=[(d,[p**j for j in range(count+1)], mult) for d,p,mult in types]
        group=[]
        for nums in splits(count,len(types)):
            d=sum(k*t[0] for k,t in zip(nums,powers))
            p=prod(t[1][k] for k,t in zip(nums,powers))
            coef=factorial(count)//prod(factorial(k) for k in nums)
            coef*=prod(t[2]**k for k,t in zip(nums,powers))
            group.append((d,p,coef))
        groups.append(group)
    for d1,p1,m1 in groups[0]:
        for d2,p2,m2 in groups[1]:
            shift=d1+d2
            p=p1*p2
            multi=m1*m2
            even=((p*free_even)&modmask).to_bytes(nbytes,'little')
            odd=((p*free_odd)&modmask).to_bytes(nbytes,'little')
            rank_e=((p*rank_even)&modmask).to_bytes(nbytes,'little')
            rank_o=((p*rank_odd)&modmask).to_bytes(nbytes,'little')
            volume=0
            for j,w in enumerate(range(shift,T+1)):
                off=(j//2)*stride
                data=odd if j&1 else even
                count=int.from_bytes(data[off:off+stride],'little')
                volume+=count
                ranks=rank_o if j&1 else rank_e
                new=int.from_bytes(ranks[off:off+stride],'little')
                cheap=((volume<<1)-count)<<a
                totals[w]+=multi*count*min(threshold,cheap,new)
    return [sum(outside[v-w]*totals[w] for w in range(v+1)) for v in range(T+1)],threshold

def certificate():
    previous = volume = success = 0
    alpha = []
    for t,bound in enumerate(full_direction_bound()):
        volume += comb(2000,t)
        k = 1+bound//Q
        current = max(previous, volume*(2*k*Q-bound)//(Q*k*(k+1)))
        assert current<=Q
        success += (current-previous)*WEIGHTS[t]
        previous = current
        alpha.append(Fraction(2,k+1))
    ml0 = 1-Fraction(success,DEN)
    return ml0, alpha

def field_checks():
    x=2
    for i in range(1,M+1):
        x=square(x)
        if i in (20,50):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert MEMORY_END<=1<<32
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


def operation_bound():
    profile_cost={}
    for width,hist in DATA:
        rows=[]
        for profile in hist:
            masses=[9**(2*width-d) for d in range(2*width+1)]
            normalizer=sum(n*w for n,w in zip(profile,masses))
            counts=Counter()
            for n,w in zip(profile,masses):
                if n:
                    counts[cost(normalizer,w)]+=n
            rows.append(tuple(sorted(counts.items())))
        profile_cost[width]=rows
    k=WORDS
    shift=3*k+3
    square=BYTES*(5+3*k)+k
    multiply=3*k+3*(shift+k)+44*k+NIBBLES*(5+shift+3*k)+k
    cube_ops=square+multiply
    anchor = 250
    leaf = 200
    widths=WIDTHS*(GROUPS-1)
    dp=[0]*(LIMIT+1)
    for i in range(len(widths)-1,-1,-1):
        width=widths[i]
        reserve=len(widths)-1-i
        inc=28
        if i==len(widths)-1:inc+=cube_ops+leaf
        elif (i+1)%len(WIDTHS)==0:inc+=cube_ops+180
        out=[0]*(LIMIT+1)
        incrow=[k+inc for k in dp]
        for row in profile_cost[width]:
            w=[0]*(LIMIT+1)
            for cc,n in row:
                for rem in range(cc+reserve,LIMIT+1):
                    w[rem]+=n*incrow[rem-cc]
            out=[max(x,y) for x,y in zip(out,w)]
        dp=[40+x for x in out]
    return 1000+sum(comb(M,dep)*(44+anchor+cube_ops+dp[LIMIT-SLOPE*dep]) for dep in range(RADIUS+1))

def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    baseline, alpha = certificate()
    run_results=[]
    for a in ((1,2),(3,4),(5,6),(7,8,9))[part]:
        nums, denom=loss(a)
        probability=sum(n*w*f for n,w,f in zip(nums,WEIGHTS,alpha))/(denom*DEN)
        run_results.append(upper_decimal(probability,10**18))
    base = upper_decimal(baseline,10**18)
    assert base <= 58309282338755975
    all_bounds = (32130851529865, 654771251756452, 1530340318922916, 1915670610504305, 1775789974200312, 1300780205703891, 705126443021130, 210429773308846, 8842704647383)
    sections = (all_bounds[:2],all_bounds[2:4],all_bounds[4:6],all_bounds[6:])
    assert all(x<=y for x,y in zip(run_results,sections[part],strict=True))
    result={"base":base,"losses":run_results}
    if part==0:
        field_checks()
        error=admission(PRECISION,LIMIT,RADIUS)
        rounded = upper_decimal(error,10**18)
        assert rounded <= 19972840061346590
        w = operation_bound()
        assert w == 70484566301335116987477129842234353111923230085460388110258
        assert w**100 <= 1<<19549
        result.update({"miss":rounded,"W":str(w),"memory_end":MEMORY_END})
    print(json.dumps(result))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
