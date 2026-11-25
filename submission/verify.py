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
    ANCHOR_BITS,ANCHOR_PARITY,MAIN_LENGTH,RADIUS,REDUCE_4,BUDGETS,BUDGET_BASE,CAP_BASE,CAP_OFFSETS,PROFILE_BASE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from dynamics import derive, data
from joint_admission import OUTPUT_SCALE, child_miss_upper
from qr import N, ORBITS, ROW_CAP, PERMUTATIONS, orbit_check, profiles as qr_compute

T=230
Q=1<<(GROUPS*M)
WEIGHTS=tuple(9**(MAIN_LENGTH-t) for t in range(T+1))
DEN=10**MAIN_LENGTH


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
    p=GROUPS*M-3
    signed=[[] for _ in range(T+1)]
    for w in range(1,p+1):
        gap=p+1-2*w
        prev,current=0,1
        for radius,row in enumerate(signed):
            row.append(abs(current))
            if radius<T:
                new,remainder=divmod(gap*current-(p-radius)*prev,radius+1)
                assert remainder==0
                prev,current=current,new
    result=prefix(MAIN_LENGTH)[1:T+2]
    for h in (0,1):
        for k in range(1-h,GROUPS):
            covered=(k+h)*M
            outside=binomial(GROUPS*M-covered)
            capacities=[]
            for w in range(1,p+1):
                value=0
                for j in range(k+1):
                    ordinary=comb(j*M,w)
                    term=comb(j*M+ANCHOR_PARITY,w)-ordinary if h else ordinary
                    value+=(-1)**(k-j)*comb(k,j)*term
                assert value>=0
                capacities.append(value)
            syst=binomial(covered)
            aux=[0]*(T+1)
            for radius,row in enumerate(signed):
                table=ranked_groups((v,c) for v,c in zip(row,capacities,strict=True) if c)
                for d in range(min(covered,T-radius)+1):
                    aux[radius+d]+=envelope(table,syst[d])
            factor=comb(GROUPS-1,k)
            for t in range(T+1):
                result[t]+=factor*sum(outside[d]*aux[t-d] for d in range(min(t,GROUPS*M-covered)+1))
    assert all(a<=b for a,b in zip(result,result[1:]))
    return [8*value for value in result]


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

def anchor_miss(alpha):
    anchor=[comb(ANCHOR_BITS,w)*9**(ANCHOR_BITS-w) if w<=ANCHOR_BITS else 0
            for w in range(T+1)]
    remain=[comb(MAIN_LENGTH-ANCHOR_BITS,w)*9**(MAIN_LENGTH-ANCHOR_BITS-w) for w in range(T+1)]
    miss=sum(anchor[d]*sum(remain[t-d]*alpha[t] for t in range(d,T+1))
             for d in range(RADIUS+1,min(ANCHOR_BITS,T)+1))
    return miss/DEN


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
def kernel(length,cut):
    parity=M*length
    ordinary=M*(length-1)
    ordinary_mult=3 if M%2==0 else 1
    volume=prefix(parity+ordinary)
    result=[volume[min(t+1,len(volume)-1)]+volume[min(t,len(volume)-1)] for t in range(cut+1)]
    rows=[[] for _ in range(cut+1)]
    for weight in range(1,parity+1):
        gap=parity+1-2*weight
        previous,current=0,1
        for radius,row in enumerate(rows):
            row.append(abs(previous+current))
            if radius<cut:
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
        coefficients=[0]*(cut+1)
        for radius,row in enumerate(rows):
            if not selected:
                coefficients[radius]+=max(row[:M])
                continue
            shifted=shifted_max(tuple(row),M)
            excluded=ranked_groups((row[w-1],count) for w,count in enumerate(capacities,1) if count)
            included=ranked_groups((shifted[w],count) for w,count in enumerate(capacities,1) if count)
            for w in range(min(inside_bits,cut-radius)+1):
                cap=comb(inside_bits,w)
                coefficients[radius+w]+=envelope(excluded,cap)+envelope(included,cap)
        choose=comb(length-1,selected)
        for t in range(cut+1):
            result[t]+=choose*sum(comb(outside_bits,b)*coefficients[t-b]
                                  for b in range(min(t,outside_bits)+1))
    return result


RUN_SCALE=10**45
CELL_BITS=76
STRIDE=(2*M+8)//8
RUN_RADICES=((83,), (81,83), (80,81), (80,), (79,80), (78,79), (78,79), (78,))
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
    for width,hist in data()+[(0,{(1,1):1})]:
        n=2*width if width else 1
        info=[]
        for h,mult in hist.items():
            denominator=100**(8*n)
            numerators=tuple((sum(q*a**(8*w)*100**(8*(n-w))
                                  for w,q in enumerate(h))*(1<<CELL_BITS)+denominator-1)
                             //denominator for a in RADICES)
            info.append((pack([mult*x for x in h],STRIDE),numerators))
        rows[width]=tuple(info)
    return rows


def class_sizes(first_pair_short):
    profiles=inputs()
    schedule=WIDTHS[:-3]+(0,)*3 if first_pair_short else WIDTHS
    widths=sorted(set(schedule),reverse=True)
    for width in widths:
        assert len(profiles[width])
    by_width=[]
    for width in widths:
        entries=profiles[width]
        count=schedule.count(width)
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
def grouped(first_pair_short):
    bins=[Counter() for _ in RADICES]
    fractional_bits=8*CELL_BITS*len(WIDTHS)
    for terms in cartesian_product(*class_sizes(first_pair_short)):
        packed=prod(block[0] for block in terms)
        for i,row in enumerate(bins):
            numerator=prod(block[1][i] for block in terms)
            k=(numerator**8-1).bit_length()-fractional_bits
            row[k]+=packed
    expected=[comb(2*M-3*first_pair_short,w) for w in range(2*M+1)]
    for rad,row in zip(RADICES,bins,strict=True):
        assert unpack(sum(row.values()),STRIDE,2*M+1)==expected
    return {rad:{k:unpack(n,STRIDE,2*M+1) for k,n in row.items()}
            for rad,row in zip(RADICES,bins,strict=True)}


def geometry_charge(length, alpha, first_pair_short, outside_short):
    p=M*length-3*outside_short
    n=(length-1)*M
    size=2*M-3*first_pair_short+p+n
    cut=min(size,T,(66,92,119,145,172,199,226,230)[length-1])
    stride=(comb(size,min(size//2,T)).bit_length()+7)//8
    # Each packed product digit is a count of a subset of all local noises.
    mask=(1<<(8*stride*(cut+1)))-1
    outside=binomial(MAIN_LENGTH-size)
    binom=binomial(p+n)
    free=pack(binom[:cut+1],stride)
    direction=kernel(length,cut)
    best=[RUN_SCALE*comb(size,w) for w in range(cut+1)]
    for a in RUN_RADICES[length-1]:
        bins=grouped(first_pair_short)[a]
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
    # Bound skipped local weights by counting every local noise.
    best.extend(RUN_SCALE*comb(size,w) for w in range(cut+1,T+1))
    channel=[a*w for a,w in zip(alpha,WEIGHTS,strict=True)]
    return sum(b*sum(outside[t-w]*channel[t] for t in range(w,T+1))
               for w,b in enumerate(best))/(RUN_SCALE*DEN)



def loss(length,alpha,shapes):
    configurations=((1,True,False),
                    (length,False,True),
                    (GROUPS-length-1,False,False))
    return sum(mult*geometry_charge(length,alpha,first,outside)
               for mult,first,outside in configurations if mult)


def baseline_and_alpha():
    previous = volume = success = 0
    alpha = []
    for t,bound in enumerate(full_direction_bound()):
        volume += comb(MAIN_LENGTH,t)
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
        if i in (3,37):
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
    assert table_word(PARITY_BASE[12]+13)==PARITY[12](13)
    assert table_word(LIST_BASE[12]+(1<<12)+5)==candidate_row(12,1)[5]
    assert table_word(LIST_BASE[1]+3)==candidate_row(1,1)[1]
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((12,0),(1,1)):
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


def main():
    part=json.load(sys.stdin)["part"]
    assert part in (0,1,2,3)
    shapes=data()
    manifest=json.loads((HERE/"submission.json").read_text())
    score=[Fraction(int(entry["score"]["numerator"]),int(entry["score"]["denominator"]))
           for entry in manifest["claims"]]
    if part==0:
        actual=qr_compute()
        assert actual==shapes[0][1]
        field_checks()
        orbits,fixed=orbit_check()
        w,_,_=derive(10)
        assert w**score[1].denominator <= 1<<score[1].numerator
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values()),
                          "W":str(w),"memory_end":MEMORY_END,"orbits":orbits,"fixed":fixed}))
        return
    declared=json.loads((HERE/"joint.json").read_text())
    joint=declared["child_miss_numerator"]
    assert type(joint) is int and 0<=joint<=OUTPUT_SCALE
    baseline,alpha=baseline_and_alpha()
    if part==1:
        assert child_miss_upper(alpha)==joint
        print(json.dumps({"child_miss_numerator":joint,"scale":OUTPUT_SCALE}))
        return
    aggregate=json.loads((HERE/"losses.json").read_text())
    assert len(aggregate)==8 and all(type(value) is int and 0<=value<=10**18 for value in aggregate)
    lengths=range(1,6) if part==2 else range(6,9)
    for length in lengths:
        measured=upper_decimal(loss(length,alpha,shapes),10**18)
        assert measured==aggregate[length-1]
    result={"runs":[[a,aggregate[a-1]] for a in lengths]}
    if part==3:
        anchor=anchor_miss(alpha)
        child=Fraction(joint,OUTPUT_SCALE)
        main=baseline+Fraction(sum(aggregate),10**18)+anchor+child
        total=1-(1-main)*Fraction(99144,100000)
        assert total<=score[0]
        result.update({"cyclic_bound":upper_decimal(main,10**18),
                       "baseline":upper_decimal(baseline,10**18),
                       "anchor":upper_decimal(anchor,10**18),
                       "child":upper_decimal(child,10**18),
                       "B_rounded":upper_decimal(total,10**18)})
    print(json.dumps(result))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
