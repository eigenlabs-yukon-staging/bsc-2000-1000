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

from pair_cube import accelerated_work, counter_samples
from dynamics import derive, data
from kernel import kernel
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
                cap=[8*sum((-1)**(k-h)*comb(k,h)*(comb(M*h+ANCHOR_BITS,w)-comb(M*h,w))
                            for h in range(k+1)) for w in range(1,1001)]
            else:
                cap=[sum((-1)**(k-h)*comb(k,h)*comb(M*h,w)
                            for h in range(k+1)) for w in range(1,1001)]
            n=M*k+ANCHOR_BITS*special
            inside=binomial(n)
            outside=binomial(1000-n)
            con=[0]*(T+1)
            for r,row in enumerate(masks):
                grouped=ranked_groups((v,c) for v,c in zip(row,cap,strict=True) if c)
                for a in range(min(n,T-r)+1):
                    con[a+r]+=envelope(grouped,(3 if special else 1)*inside[a])
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

DEFECT=M-ANCHOR_BITS
RUN_SCALE=1<<224
BITS=76
RADIX=80
BINS=16
PAIR_CUT=70
PAIR_STRIDE=(comb(2*M,PAIR_CUT).bit_length()+7)//8
PAIR_MASK=(1<<(8*PAIR_STRIDE*(PAIR_CUT+1)))-1

def pack(row,stride):
    return int.from_bytes(b''.join(x.to_bytes(stride,'little') for x in row),'little')

def unpack(value,stride,count):
    buf=value.to_bytes(stride*count,'little')
    return [int.from_bytes(buf[i:i+stride],'little') for i in range(0,len(buf),stride)]

def step(register,bits,hist,rad):
    options=Counter()
    for h,mult in hist.items():
        denominator=100**(8*bits)
        numerator=sum(n*rad**(8*w)*100**(8*(bits-w)) for w,n in enumerate(h))
        top=((numerator<<BITS)+denominator-1)//denominator
        k=(top**BINS-1).bit_length()-BITS*BINS
        options[k]+=pack([mult*n for n in h],PAIR_STRIDE)
    out=Counter()
    for k,x in register.items():
        for d,y in options.items():
            out[k+d]+=x*y
    return {k:n&PAIR_MASK for k,n in out.items() if n&PAIR_MASK}

@lru_cache(None)
def common(rad):
    state={0:1}
    hist=dict(data())
    for width in WIDTHS[:-DEFECT]:
        state=step(state,2*width,hist[width],rad)
    return state

@lru_cache(None)
def group(kind,rad):
    state=common(rad)
    if kind==0:
        bits1,hist1=1,{(1,0):1,(0,1):1}
    elif kind==1:
        bits1,hist1=1,{(1,1):1}
    else:
        bits1,hist1=2,dict(data())[1]
    for _ in range(DEFECT):
        state=step(state,bits1,hist1,rad)
    n=2*M-DEFECT*(kind<2)
    assert unpack(sum(state.values()),PAIR_STRIDE,PAIR_CUT+1)==[comb(n,i) for i in range(PAIR_CUT+1)]
    return {k:unpack(p,PAIR_STRIDE,PAIR_CUT+1) for k,p in state.items()}

@lru_cache(None)
def root_up(remainder):
    n=1<<(1024+remainder)
    f=isqrt(isqrt(isqrt(isqrt(n))))
    return f+(f**16<n)

def proper_geometry_counts(length,kind,p,free,local,has_special,rad):
    cut=(84,124,161,198,230,230)[length-1]
    theory=kernel(M,DEFECT,length,kind==0,has_special and kind!=0,cut)
    stride=(comb(local,min(cut,local//2)).bit_length()+7)//8
    free_packed=pack(binomial(free)[:cut+1],stride)
    mask=(1<<(8*stride*(cut+1)))-1
    bounds=[]
    max_n=0
    rp=1
    for d in theory[:cut+1]:
        max_n=max(d*rp,max_n*100**8)
        den=rp<<(p+1)
        bounds.append((max_n*RUN_SCALE+den-1)//den)
        rp*=rad**8
    by_weight=[0]*(cut+1)
    for k,h in group(kind,rad).items():
        q,j=divmod(k,BINS)
        top=root_up(j)
        shift=64-q
        add=(1<<shift)-1 if shift>=0 else 0
        nums=unpack((pack(h[:cut+1],stride)*free_packed)&mask,stride,cut+1)
        for w,n in enumerate(nums):
            if n:
                mass=bounds[w]*top
                prob=(mass+add)>>shift if shift>=0 else mass<<-shift
                by_weight[w]+=n*min(RUN_SCALE,prob)
    observed=2*M-DEFECT*(kind<2)
    freeout=binomial(2000-local)
    for w in range(cut+1):
        extra=sum(comb(observed,d)*comb(free,w-d)
                  for d in range(PAIR_CUT+1,min(observed,w)+1))
        by_weight[w]=min(RUN_SCALE*comb(local,w),by_weight[w]+RUN_SCALE*extra)
    by_weight.extend(RUN_SCALE*comb(local,w) for w in range(cut+1,T+1))
    return [Fraction(sum(by_weight[w]*freeout[t-w] for w in range(t+1)),RUN_SCALE)
            for t in range(T+1)]


def geometry(start,length):
    fields=[(start+i)%GROUPS for i in range(length+1)]
    parity=sum(PAR_WIDTHS[i] for i in fields[1:])
    free=parity+sum(SYS_WIDTHS[i] for i in fields[1:-1])
    local=free+SYS_WIDTHS[start]+PAR_WIDTHS[start]
    return min(start,2),parity,free,local,0 in fields[:-1]


def scaled_proper_counts(lengths):
    totals=[Fraction(0)]*(T+1)
    for a in lengths:
        shapes=Counter(geometry(start,a) for start in range(GROUPS))
        for key,mult in shapes.items():
            row=proper_geometry_counts(a,*key,RADIX)
            totals=[x+mult*y for x,y in zip(totals,row,strict=True)]
    numerators=[value*RUN_SCALE for value in totals]
    assert all(value.denominator==1 for value in numerators)
    return [value.numerator for value in numerators]


def block_error(child,proper):
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
    for j in range(1,M+1):
        x=square(x)
        if j in (11,13):
            a,b=x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a==1
    assert x==2
    assert sum(SYS_WIDTHS)==sum(PAR_WIDTHS)==1000
    assert len(SEARCH_WIDTHS)==TOTAL_CELLS==66
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
        _,caps,_=derive(10,True)
        w = accelerated_work(caps)
        assert w == 9696361624460279422229763503205612703275629028686308275637864
        assert w**6695 <= 1<<1356361
        coverage = counter_samples()
        assert coverage["cube_ops"] == 712
        cutoffs=[217,218,219]
        pairs=pair_rows()
        for t in cutoffs:
            assert pair_direction_bound(t)==pairs[t]
        print(json.dumps({"syndrome_types":len(actual),"syndromes":sum(actual.values()),
                          "checked_cutoffs":cutoffs,
                          "W":str(w),"orbits":orbits,"fixed":fixed,**coverage}))
        return
    if part==1:
        rows=pair_rows()
        cutoffs=list(range(220,227))
        for t in cutoffs:
            assert pair_direction_bound(t)==rows[t]
        print(json.dumps({"checked_cutoffs":cutoffs}))
        return
    path=json.loads((HERE/"misses.json").read_text())
    assert isinstance(path,list) and len(path)==T+1
    child=[int(value) for value in path]
    assert all(0<=value<=comb(2000,t) for t,value in enumerate(child))
    rows=json.loads((HERE/"proper.json").read_text())
    assert isinstance(rows,list) and len(rows)==T+1
    early=[int(value) for value in rows]
    assert all(0<=value<=(GROUPS*2)*RUN_SCALE*comb(2000,t)
               for t,value in enumerate(early))
    if part==2:
        assert child_removed_counts(T)==child
        assert scaled_proper_counts(range(1,3))==early
        print(json.dumps({"checked_weights":len(child),"checked_lengths":list(range(1,3))}))
    else:
        late=scaled_proper_counts(range(3,GROUPS))
        proper=[(a+b+RUN_SCALE-1)//RUN_SCALE for a,b in zip(early,late,strict=True)]
        error=block_error(child,proper)
        assert error<=Fraction(5808462030158777248264341175977,100000000000000000000000000000000)
        print(json.dumps({"B_rounded":upper_decimal(error,10**32),"scale":10**32}))


def upper_decimal(value,scale):
    return -(-value.numerator*scale//value.denominator)


if __name__=='__main__':
    main()
