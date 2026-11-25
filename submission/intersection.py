"""Exact signed parity-shell bounds for full-change Hamming intersections."""
from bisect import bisect_left
from functools import lru_cache
from math import comb

from construction import ANCHOR_BITS, GROUPS, M

K=1000
N=2000
DEFECT=M-ANCHOR_BITS

def pack(values,stride):
    return int.from_bytes(b''.join(value.to_bytes(stride,'little') for value in values),'little')


def unpack(value,stride,count):
    raw=value.to_bytes(stride*count,'little')
    return [int.from_bytes(raw[i:i+stride],'little') for i in range(0,len(raw),stride)]


def overlap(t,d):
    full=t-(d+1)//2
    if full<0:
        return 0
    tail=[]
    total=0
    for j in range(full+1):
        total+=comb(N-d,j)
        tail.append(total)
    return sum(comb(d,a)*tail[t-max(a,d-a)]
               for a in range(max(0,d-t),min(d,t)+1))


def rank(magnitudes,capacities):
    values=sorted(((x,c) for x,c in zip(magnitudes,capacities,strict=True) if c),reverse=True)
    counts=[0]
    total=[0]
    for value,count in values:
        counts.append(counts[-1]+count)
        total.append(total[-1]+value*count)
    return values,counts,total


def top(table,n):
    rows,counts,total=table
    if n>=counts[-1]:
        return total[-1]
    i=bisect_left(counts,n)-1
    return total[i]+(n-counts[i])*rows[i][0]


@lru_cache(None)
def shape_counts():
    ordinary_mu=1 if M%2 else 3
    shortened_mu=4<<DEFECT
    ordinary=[[comb(j*M,d) for d in range(K+1)] for j in range(GROUPS)]
    special=[[comb(j*M+ANCHOR_BITS,d)-ordinary[j][d] for d in range(K+1)]
             for j in range(GROUPS)]
    counts={}
    bins={}
    for s in (0,1):
        for k in range(GROUPS):
            sign=[(-1)**(k-j)*comb(k,j) for j in range(k+1)]
            rows=special if s else ordinary
            row=[sum(sign[j]*rows[j][d] for j in range(k+1)) for d in range(K+1)]
            assert all(value>=0 for value in row)
            counts[k,s]=row
            bins[k,s]=[ordinary_mu**k*shortened_mu**s*n for n in row[1:]]
    return counts,bins


def exact_columns(t):
    degree=2*t
    intersections=[overlap(t,d) for d in range(degree+1)]
    V=intersections[0]
    stride=(V.bit_length()+K+8)//8
    mask=(1<<(8*stride*(degree+1)))-1
    h=pack(reversed(intersections),stride)
    shell=pack((comb(K,p) for p in range(degree+1)),stride)
    row=unpack((h*shell)&mask,stride,degree+1)
    columns=[]
    for _ in range(K):
        previous=row
        row=[previous[0]]
        for p in range(1,degree+1):
            row.append(previous[p]-previous[p-1]-row[-1])
        columns.append([abs(n) for n in row])
    return V,columns


def pair_direction_bound(t):
    counts,bins=shape_counts()
    V,columns=exact_columns(t)
    extra=0
    eta=(1<<(DEFECT+1))-1
    for d in range(GROUPS,2*t+1):
        magnitudes=[column[2*t-d] for column in columns]
        for s in (0,1):
            for k in range(1-s,GROUPS):
                table=rank(magnitudes,bins[k,s])
                outside_k=GROUPS-1-k
                full_outside=M*outside_k+ANCHOR_BITS*(1-s)
                q_min=max(outside_k+1-s,d-(k*M+s*ANCHOR_BITS))
                q_max=min(full_outside,d-k-s)
                extra+=comb(GROUPS-1,k)*sum(
                    counts[outside_k,1-s][q]*top(table,eta**s*counts[k,s][d-q])
                    for q in range(q_min,q_max+1)
                )
    return (V*V+extra+V-1)//V
