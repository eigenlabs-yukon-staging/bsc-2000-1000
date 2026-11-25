"""Exact full-ML, proper-run, and restricted-search certificates."""
from bisect import bisect_left
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import json
from itertools import combinations_with_replacement, product as cartesian_product
from math import comb, isqrt, prod
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
def kernel(length):
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


RUN_SCALE=10**45
CELL_BITS=76
STRIDE=(2*M+8)//8
RUN_RADICES=((83,), (81,83), (80,81), (80,), (79,80), (78,79), (78,79))
RADICES=tuple(sorted({r for row in RUN_RADICES for r in row}))


def pack(row,stride):
    return int.from_bytes(b''.join(x.to_bytes(stride,'little') for x in row),'little')


def unpack(value,stride,count):
    packed=value.to_bytes(stride*count,'little')
    return [int.from_bytes(packed[i:i+stride],'little') for i in range(0,len(packed),stride)]


def upper_pow(frac_n,frac_d):
    power=frac_n<<512
    r=isqrt(isqrt(isqrt(power//frac_d)))
    return r+(r**8*frac_d<power)


@lru_cache(None)
def inputs():
    rows={}
    for width,hist in data():
        n=2*width
        info=[]
        for h,mult in hist.items():
            denominator=100**(8*n)
            numerators=tuple((sum(q*a**(8*w)*100**(8*(n-w))
                                  for w,q in enumerate(h))*(1<<CELL_BITS)+denominator-1)
                             //denominator for a in RADICES)
            info.append((pack([mult*x for x in h],STRIDE),numerators))
        rows[width]=tuple(info)
    return rows


def class_sizes():
    profiles=inputs()
    widths=sorted(set(WIDTHS),reverse=True)
    for width in widths:
        assert len(profiles[width])
    by_width=[]
    for width in widths:
        entries=profiles[width]
        count=WIDTHS.count(width)
        rows=[]
        for selected in combinations_with_replacement(range(len(entries)),count):
            factor=1
            for j in range(2,count+1):
                factor*=j
            for repeat in Counter(selected).values():
                for j in range(2,repeat+1):
                    factor//=j
            packed=prod(entries[index][0] for index in selected)*factor
            numerators=tuple(prod(entries[index][1][j] for index in selected)
                             for j in range(len(RADICES)))
            rows.append((packed,numerators))
        by_width.append(rows)
    return by_width


@lru_cache(None)
def grouped():
    bins=[Counter() for _ in RADICES]
    fractional_bits=8*CELL_BITS*len(WIDTHS)
    for terms in cartesian_product(*class_sizes()):
        packed=prod(block[0] for block in terms)
        for i,row in enumerate(bins):
            numerator=prod(block[1][i] for block in terms)
            k=(numerator**8-1).bit_length()-fractional_bits
            row[k]+=packed
    expected=[comb(2*M,w) for w in range(2*M+1)]
    for rad,row in zip(RADICES,bins,strict=True):
        assert unpack(sum(row.values()),STRIDE,2*M+1)==expected
    return {rad:{k:unpack(n,STRIDE,2*M+1) for k,n in row.items()}
            for rad,row in zip(RADICES,bins,strict=True)}


def loss(length, alpha, data):
    p=M*length
    n=(length-1)*M
    size=2*M+p+n
    cut=min(size,T)
    stride=(comb(size,min(size//2,T)).bit_length()+7)//8
    # Each packed product digit is a count of a subset of all local noises.
    mask=(1<<(8*stride*(cut+1)))-1
    outside=binomial(2000-size)
    binom=binomial(p+n)
    free=pack(binom[:cut+1],stride)
    direction=kernel(length)
    best=[RUN_SCALE*comb(size,w) for w in range(cut+1)]
    for a in RUN_RADICES[length-1]:
        bins=grouped()[a]
        max_num=0
        ap=1
        factor=[]
        for d in direction[:cut+1]:
            max_num=max(d*ap,max_num*100**8)
            factor.append((max_num*RUN_SCALE,ap<<(p+1)))
            ap*=a**8
        loss=[0]*(cut+1)
        for k,hist in bins.items():
            exponent,root=divmod(k,8)
            up=upper_pow(1<<root,1)
            if exponent>=0:
                scale_n=up<<exponent
                scale_d=1<<64
            else:
                scale_n=up
                scale_d=1<<(64-exponent)
            mass=unpack((pack(hist[:cut+1],stride)*free)&mask,stride,cut+1)
            for w,frequency in enumerate(mass):
                if frequency:
                    num,den=factor[w]
                    den*=scale_d
                    prob=min(RUN_SCALE,(num*scale_n+den-1)//den)
                    loss[w]+=prob*frequency
        best=[min(x,y) for x,y in zip(best,loss,strict=True)]
    channel=[a*w for a,w in zip(alpha,WEIGHTS,strict=True)]
    return sum(b*sum(outside[t-w]*channel[t] for t in range(w,T+1))
               for w,b in enumerate(best))/(RUN_SCALE*DEN)



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
        assert w == 1891535973990986868106847676987599045627015342793609452
        assert w**125 <= 1<<22538
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
        assert total <= Fraction(80289546035962436,10**18)
        print(json.dumps({"B_rounded":upper_decimal(total,10**18),
            "baseline_rounded":upper_decimal(baseline,10**18),
            "run_losses_rounded":[upper_decimal(x,10**18) for x in losses],
            "anchor_charge_rounded":upper_decimal(anchor,10**18),
            "child_charge_rounded":upper_decimal(child,10**18)}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
