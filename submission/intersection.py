"""Integer full-change Hamming-ball intersection bound for two degree-500 fields."""
from bisect import bisect_left
from math import comb


M=500
N=2000
K=1000


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

def pair_direction_bound(t):
    degree=2*t
    intersections=[overlap(t,d) for d in range(degree+1)]
    volume=intersections[0]
    stride=(volume.bit_length()+K+8)//8
    mask=(1<<(8*stride*(degree+1)))-1
    left=pack(reversed(intersections),stride)
    first=[0]*len(MASK_BINS)
    full=[0]*len(MASK_BINS)
    for w in range(1,K+1):
        j=bisect_left(MASK_BINS,min(w,K-w)+1)-1
        first[j]+=3*comb(M,w)
        full[j]+=9*(comb(K,w)-2*comb(M,w))
    columns=[]
    for j in MASK_BINS:
        even=pack((comb(j,p//2) if p%2==0 else 0
                   for p in range(min(degree,2*j)+1)),stride)
        linear=pack((comb(K-2*j,p) for p in range(min(degree,K-2*j)+1)),stride)
        packed=left*((even*linear)&mask)&mask
        columns.append(unpack(packed,stride,degree+1))
    correction=0
    cM=[comb(M,d) for d in range(degree+1)]
    for d in range(2,degree+1):
        mags=[column[degree-d] for column in columns]
        assert all(a>=b for a,b in zip(mags,mags[1:]))
        rank_one=rank(mags,first)
        correction+=2*sum(cM[b]*top(rank_one,cM[d-b]) for b in range(1,d))
        correction+=top(rank(mags,full),comb(K,d)-2*cM[d])
    D=(volume*volume+correction+volume-1)//volume
    return D
