from bisect import bisect_left
from collections import deque
from functools import lru_cache
from math import comb

T=230

def grouped(rows):
    ordered=sorted(rows,reverse=True)
    counts=[0];totals=[0]
    for value,mult in ordered:
        counts.append(counts[-1]+mult)
        totals.append(totals[-1]+mult*value)
    return ordered,counts,totals

def top(table,n):
    rows,counts,totals=table
    n=min(n,counts[-1])
    j=bisect_left(counts,n)-1
    return totals[j]+(n-counts[j])*rows[j][0]

def shifted(row,width):
    q=deque()
    out=[]
    for i,x in enumerate(row):
        while q and q[-1][1]<=x:q.pop()
        q.append((i,x))
        if i>=width-1:
            while q[0][0]<=i-width:q.popleft()
            out.append(q[0][1])
    return out

@lru_cache(None)
def signs(p):
    rows=[[] for _ in range(T+1)]
    for w in range(1,p+1):
        g=p+1-2*w
        prev,actual=0,1
        for l,row in enumerate(rows):
            row.append(abs(actual))
            if l<T:
                new,rest=divmod(g*actual-(p-l)*prev,l+1)
                assert rest==0
                prev,actual=actual,new
    return rows

def kernel(m,defect,a,first_short,other_short):
    assert not(first_short and other_short)
    short=m-defect
    first_width=short if first_short else m
    first_count=(1<<(defect+1))-1 if first_short else 1
    other_count=a-1-other_short
    ordinary_mult=1 if m%2 else 3
    p=a*m-defect*(first_short+other_short)
    n=(a-1)*m-defect*other_short
    volume=[];total=0
    for t in range(T+1):
        total+=comb(n+p,t);volume.append(total)
    for h in range(1+other_short):
        for k in range(other_count+1):
            degree=k*m+h*short
            outside=tuple(comb(n-degree,d) for d in range(T+1))
            formal=tuple(comb(degree,d) for d in range(T+1))
            capacity=[]
            mask_mult=ordinary_mult**k*((4<<defect) if h else 1)
            direction_mult=((1<<(defect+1))-1 if h else 1)
            for w in range(1,degree+1):
                cap=0
                for l in range(k+1):
                    val=comb(m*l,w)
                    if h:
                        val=comb(m*l+short,w)-val
                    cap+=(-1)**(k-l)*comb(k,l)*val
                capacity.append(cap*mask_mult)
            aux=[0]*(T+1)
            for r,row in enumerate(signs(p)):
                if degree==0:
                    aux[r]+=first_count*max(row[:first_width])
                    continue
                slide=shifted(row,first_width)
                listed=grouped((row[w-1],num) for w,num in enumerate(capacity,1) if num)
                joined=grouped((slide[w],num) for w,num in enumerate(capacity,1) if num)
                for d in range(min(degree,T-r)+1):
                    choices=direction_mult*formal[d]
                    aux[r+d]+=top(listed,choices)+first_count*top(joined,choices)
            mult=comb(other_count,k)
            for t in range(T+1):
                volume[t]+=mult*sum(outside[r]*aux[t-r] for r in range(min(n-degree,t)+1))
    return volume
