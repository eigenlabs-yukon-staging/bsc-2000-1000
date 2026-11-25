"""Integer full-change Hamming-ball intersection bounds."""
from bisect import bisect_left
from functools import lru_cache
from math import comb

from construction import GROUPS,M

K=GROUPS*M
N=2*K

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


MASK_BINS=(0,1,2,3,4,6,8,12,16,24,32,48,64,80,96,112,128,144,160,176,192,208,224,
           240,256,280,304,328,352,376,400,424,448,472,500)

@lru_cache(None)
def shape_counts():
    g=GROUPS
    mu=1 if M%2 else 3
    mask_bins=tuple(x for x in MASK_BINS if x<=K//2)
    all_binom=[[comb(j*M,d) for d in range(K+1)] for j in range(g+1)]
    multiplicities=[[] for _ in range(g+1)]
    capacities=[[] for _ in range(g+1)]
    for k in range(g+1):
        weights=[(-1)**(k-j)*comb(k,j) for j in range(k+1)]
        multiplicities[k]=[sum(weights[j]*all_binom[j][w] for j in range(k+1))
                           for w in range(K+1)]
        assert all(n>=0 for n in multiplicities[k])
        capacities[k]=[0]*len(mask_bins)
        for w in range(1,K+1):
            bucket=bisect_left(mask_bins,min(w,K-w)+1)-1
            capacities[k][bucket]+=mu**k*multiplicities[k][w]
    return mask_bins,multiplicities,capacities


def pair_direction_bound(t):
    g=GROUPS
    degree=2*t
    mask_bins,multiplicities,capacities=shape_counts()
    J=[overlap(t,d) for d in range(degree+1)]
    V=J[0]
    stride=(V.bit_length()+K+8)//8
    mask=(1<<(8*stride*(degree+1)))-1
    packed_J=pack(reversed(J),stride)
    columns=[]
    for j in mask_bins:
        u=pack((comb(j,p//2) if p%2==0 else 0 for p in range(min(degree,2*j)+1)),stride)
        h=pack((comb(K-2*j,p) for p in range(min(degree,K-2*j)+1)),stride)
        columns.append(unpack((packed_J*((u*h)&mask))&mask,stride,degree+1))

    correction=0
    for d in range(g,degree+1):
        magnitudes=[col[degree-d] for col in columns]
        for k in range(1,g+1):
            limiting=rank(magnitudes,capacities[k])
            b0=max(g-k,d-k*M)
            b1=min((g-k)*M,d-k)
            correction+=comb(g,k)*sum(multiplicities[g-k][b]*
                                       top(limiting,multiplicities[k][d-b])
                                       for b in range(b0,b1+1))
    return (V*V+correction+V-1)//V

