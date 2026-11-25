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
    ANCHOR_BITS,RADIUS,REDUCE_4,BUDGETS,BUDGET_BASE,CAP_BASE,CAP_OFFSETS,PROFILE_BASE,TOTAL_CELLS,candidate_row,cost,encode,
    poly_remainder,square,table_word,BYTES,NIBBLES,WIDTHS,WORDS,QR_ROW_BASE,
    SQUARE_BASE,PARITY_BASE,LIST_BASE,QR_PARITY_BASE)

from pair_cube import accelerated_work, counter_samples
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
RUN_RADICES=((82,83,84),(80,81,82,83),(79,80,81,82),(79,80,81),(78,79,80,81),(78,79,80),(77,78,79,80),(77,78,79),(77,78,79))
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


def proper_length_counts(length):
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
    return [Fraction(sum(best[w]*outside[t-w] for w in range(min(t,cut)+1)),RUN_SCALE)
            for t in range(T+1)]


def scaled_proper_counts(lengths):
    rows=[proper_length_counts(a) for a in lengths]
    scaled=[GROUPS * RUN_SCALE * sum(items) for items in zip(*rows, strict=True)]
    assert all(value.denominator == 1 for value in scaled)
    return [value.numerator for value in scaled]


def block_error(child, proper):
    anchor=anchor_counts()
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
        if i in (20,50):
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
    assert table_word(QR_PARITY_BASE+734)==PARITY[24](734)
    assert max(table_word(REDUCE_4+i) for i in range(16)) < 1<<64
    for width,syndrome in ((4,s) for s in range(16)):
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
        field_checks()
        orbits,fixed=orbit_check()
        _,caps,_=derive(10,True)
        w=accelerated_work(caps)
        assert w == 21306706949390502987836000332430252282142597733
        assert w**10 <= 1<<1539
        coverage=counter_samples()
        assert coverage["cube_ops"]==227
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values()),
                          "W":str(w),"orbits":orbits,"fixed":fixed,**coverage}))
        return
    early=list(map(int,json.loads((HERE/"proper.json").read_text())))
    assert len(early)==T+1
    assert all(0<=count<=5*GROUPS*RUN_SCALE*comb(2000,t)
               for t,count in enumerate(early))
    if part==1:
        rows=pair_rows()
        for t,declared in rows.items():
            assert pair_direction_bound(t)==declared
        assert scaled_proper_counts(range(1,6))==early
        print(json.dumps({"checked_cutoffs":list(rows),"checked_proper_lengths":[1,2,3,4,5]}))
        return
    path=json.loads((HERE/"misses.json").read_text())
    assert isinstance(path,list) and len(path)==T+1
    child=[int(value) for value in path]
    assert all(0<=value<=comb(2000,t) for t,value in enumerate(child))
    if part==2:
        assert child_removed_counts(T)==child
        print(json.dumps({"checked_weights":len(child)}))
    else:
        late=scaled_proper_counts(range(6,GROUPS))
        proper=[(a+b+RUN_SCALE-1)//RUN_SCALE for a,b in zip(early,late,strict=True)]
        error=block_error(child,proper)
        assert error<=Fraction(1404719304643696798689007504653,10000000000000000000000000000000)
        print(json.dumps({"B_rounded":upper_decimal(error,10**32),"scale":10**32,
                          "checked_proper_lengths":list(range(6,GROUPS))}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
