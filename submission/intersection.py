"""Exact signed parity-shell bounds for full-change Hamming intersections."""
from bisect import bisect_left
from functools import lru_cache
from math import comb

from construction import GROUPS, M

G=GROUPS
K=G*M
N=2*K
mu=1 if M%2 else 3

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


def rank(magnitudes, capacities, order):
    values=[(magnitudes[i],capacities[i]) for i in order if capacities[i]]
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
def counts():
    binom=[[comb(j*M,w) for w in range(K+1)] for j in range(G+1)]
    multiplicities=[]
    for k in range(G+1):
        signs=[(-1)**(k-j)*comb(k,j) for j in range(k+1)]
        row=[sum(signs[j]*binom[j][w] for j in range(k+1)) for w in range(K+1)]
        assert all(n>=0 for n in row)
        multiplicities.append(row)
    return multiplicities, [[mu**k*n for n in row[1:]]
                            for k,row in enumerate(multiplicities)]


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
    multiplicities,cap=counts()
    V,columns=exact_columns(t)
    correction=0
    for d in range(G,2*t+1):
        magnitudes=[column[2*t-d] for column in columns]
        order=sorted(range(len(magnitudes)),key=magnitudes.__getitem__,reverse=True)
        for k in range(1,G+1):
            table=rank(magnitudes,cap[k],order)
            complement=multiplicities[G-k]
            chosen=multiplicities[k]
            q0=max(G-k,d-k*M)
            q1=min((G-k)*M,d-k)
            correction+=comb(G,k)*sum(complement[q]*top(table,chosen[d-q])
                                       for q in range(q0,q1+1))
    return (V*V+correction+V-1)//V


