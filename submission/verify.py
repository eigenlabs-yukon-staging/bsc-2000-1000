"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter, deque
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
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from dynamics import derive, data, fixed_tails
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
            sum((-1)**(k-h)*comb(k,h)*comb(h*M,w) for h in range(1,k+1))
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

SCALE=10**60
NOISE_DENOMINATORS=(10,12)


def weighted_miss(tails,alpha):
    anchor=[comb(M,w)*9**(M-w) if w<=M else 0 for w in range(T+1)]
    outside=anchor
    remain=[comb(2000-M,w)*9**(2000-M-w) for w in range(T+1)]
    anchor_fail=Fraction(0)
    for d in range(RADIUS+1,min(M,T)+1):
        anchor_fail+=anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
    anchor_fail/=DEN
    diffs=[x-y for x,y in zip(alpha,alpha[1:]+[Fraction(0)],strict=True)]
    g=[]
    for s in range(T+1):
        p=sum(outside[w]*diffs[s+w] for w in range(T-s+1))/10**M
        g.append(upper_decimal(p,SCALE))
    ratios={}
    for den in NOISE_DENOMINATORS[1:]:
        ratios[den]=[(den**1750*9**(1750-z),10**1750*(den-1)**(1750-z))
                     for z in range(T+1)]
    child_fail=0
    for d in range(RADIUS+1):
        subtotal=0
        for z in range(T-d+1):
            best=tails["10"][d]
            for den in NOISE_DENOMINATORS[1:]:
                numerator,denom=ratios[den][z]
                adjusted=(tails[str(den)][d]*numerator+denom-1)//denom
                best=min(best,adjusted)
            subtotal+=g[d+z]*best
        child_fail+=anchor[d]*subtotal
    return anchor_fail,Fraction(child_fail,10**M*SCALE**2)


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
        multi=[sum((-1)**(k-h)*comb(k,h)*comb(h*m,w)
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
                total[l+x]+=envelope(entries_nonfirst,cap)+envelope(entries_withfirst,cap)
        choose=comb(a-1,k)
        for n in range(T+1):
            result[n]+=choose*sum(outside[b]*total[n-b] for b in range(min(n,ordinary-degree)+1))
    return result


GRID=(((83, 1), (83, 8), (83, 4), (83, 6)), ((81, 3), (80, 1), (80, 2), (81, 5), (81, 4), (81, 1), (81, 6)), ((80, 3), (79, 1), (80, 1), (80, 2), (80, 4), (80, 5)), ((79, 2), (79, 1), (79, 3), (79, 4), (79, 5)), ((79, 3), (79, 2), (79, 1), (78, 1), (79, 4), (79, 5)), ((78, 2), (78, 3), (78, 1), (79, 4), (78, 4), (79, 5)), ((78, 4), (78, 3), (78, 5), (78, 2), (78, 6)))


def power_up(n,d,j):
    numerator=n**j<<512
    denominator=d**j
    floor=isqrt(isqrt(isqrt(numerator//denominator)))
    return floor+(floor**8*denominator<numerator)


def run_coefficients(length, a, b, j, data):
    free=M*(2*length-1)
    parts=[]
    for (width,hist),count in zip(data,(5,1,1),strict=True):
        row=[0]*(2*width+1)
        for profile,multiplicity in hist.items():
            n=sum(v*a**(8*i)*b**(16*width-8*i) for i,v in enumerate(profile))
            root=power_up(n,b**(16*width),j)
            for i,v in enumerate(profile):
                row[i]+=multiplicity*v*root
        parts.append((row,count))
    parts.append(([comb(free,i) for i in range(T+1)],1))
    row=product(parts,T,4000)
    exponent,rem=divmod(j*free,8)
    top,bottom=a**8+b**8,b**8
    e,rem2=divmod(j*(M-1)*length,8)
    numerator=top**exponent*power_up(top,bottom,rem)*power_up(1,1<<rem2,1)
    denominator=bottom**exponent*(1<<(e+(7+2)*64))
    plain=Fraction(numerator,denominator)
    e,rem2=divmod(j*M*length,8)
    scale=Fraction(power_up(1,1<<rem2,1),1<<(e+(7+2)*64))
    cumulative=Fraction(0)
    projected=Fraction(0)
    out=[]
    direction=direction_bound(M,length)
    for w,n in enumerate(row):
        current=Fraction(direction[w]*a**(8*w),b**(8*w))
        if current>cumulative:
            cumulative=current
            projected=scale*power_up(cumulative.numerator,cumulative.denominator,j)
        out.append(min(plain,projected)*Fraction(n*b**(j*w),a**(j*w)))
    return out


def loss(length, alpha, data):
    local = (2*length+1)*M
    by_weight = [Fraction(comb(local,w)) for w in range(T+1)]
    for a,j in GRID[length-1]:
        by_weight = [min(x,y) for x,y in zip(by_weight,run_coefficients(length,a,100,j,data),strict=True)]
    outside = binomial(2000-local)
    tilted = [w*g for w,g in zip(WEIGHTS,alpha,strict=True)]
    propagate = [sum(outside[t-w]*tilted[t] for t in range(w,T+1)) for w in range(T+1)]
    return sum(x*y for x,y in zip(by_weight,propagate,strict=True))/DEN


def certificate(data):
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
    losses = [loss(length,alpha,data) for length in range(1,GROUPS)]
    return ml0, losses, alpha


def field_checks():
    x=2
    for i in range(1,M+1):
        x=square(x)
        if i==25:
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert QR_ROW_BASE + ORBITS*ROW_CAP == BUDGET_BASE
    assert BUDGET_BASE+len(BUDGETS)==PROFILE_BASE[24]
    assert CAP_OFFSETS[0]==CAP_BASE and CAP_OFFSETS[-1]==MEMORY_END<=1<<32
    assert len(BUDGETS)==RADIUS+1 and min(BUDGETS)>=TOTAL_CELLS
    assert all(table_word(BUDGET_BASE+d)==c for d,c in enumerate(BUDGETS))
    assert table_word(SQUARE_BASE)==0
    assert table_word(SQUARE_BASE+3*WORDS)==square(3)&((1<<64)-1)
    assert table_word(PARITY_BASE[4]+13)==PARITY[4](13)
    assert table_word(LIST_BASE[4]+21)==candidate_row(4,1)[5]
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((4,0),(4,15)):
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
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values())}))
        return
    tails=json.loads((HERE/"tails.json").read_text())
    assert set(tails)=={str(q) for q in NOISE_DENOMINATORS}
    assert all(len(row)==RADIUS+1 and all(0<=x<=SCALE for x in row) for row in tails.values())
    if part==1:
        field_checks()
        orbits,fixed=orbit_check()
        w,_,nominal=derive(10)
        assert w == 538712329465928215968257672196566561182176754236
        assert w**1000 <= 1<<158561
        assert fixed_tails(nominal)==tails["10"]
        print(json.dumps({"W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
    elif part==2:
        for den in NOISE_DENOMINATORS[1:]:
            _,_,success=derive(den)
            assert fixed_tails(success)==tails[str(den)]
        print(json.dumps({"checked_tilts":list(NOISE_DENOMINATORS[1:])}))
    else:
        baseline,losses,alpha=certificate(shapes)
        anchor,child=weighted_miss(tails,alpha)
        total=baseline+GROUPS*sum(losses)+anchor+child
        assert total <= Fraction(193785795755191728,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
