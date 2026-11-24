"""Verify every selected minor at a fixed witness for one orientation."""
from itertools import combinations
from pathlib import Path
import json
import sys

HERE=Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from construction import WIDTH, chart_matrices, remainder

WITNESSES=(0xae03b85f10141192e441d5899bcba22ffc91ad0f845cd8f326f282c9b38d00338e7c416a39a3ba8e8b08a65c7fdb6619f5fa301c6ac0a73f7fa9b3030acfb,0x1d76845fc782cd4f80b7b1bbbfcb02124bff11c2253850054d49f3fc3544d1c367f0c85d9cd72139b7911801efb83fb0e0c0e20197775ae2761492a07e0a5)

def build_lut(v):
    columns=[v]
    for _ in range(1,WIDTH):
        columns.append(remainder(columns[-1]<<1))
    values=[]
    for i in range(7):
        t=[0]
        for j in range(8):
            p=8*i+j
            x=columns[p] if p<WIDTH else 0
            t.extend(v^x for v in t[:])
        values.append(tuple(t))
    return tuple(values)

def times(t,x):
    return (t[0][x&255]^t[1][x>>8&255]^t[2][x>>16&255]^t[3][x>>24&255]
            ^t[4][x>>32&255]^t[5][x>>40&255]^t[6][x>>48])

def compute(matrix):
    indexed=[[(sum(1<<x for x in c),c) for c in combinations(range(10),k)]
             for k in range(7)]
    expansion=[tuple((mask,tuple((j,mask^(1<<j)) for j in cols))
                     for mask,cols in stage) for stage in indexed]
    luts={v:build_lut(v) for v in set(x for row in matrix for x in row)}
    matrices=tuple(tuple(luts[x] for x in row) for row in matrix)

    def extend(previous,r,k):
        row=matrices[r]
        output=[0]*1024
        for mask,terms in expansion[k]:
            v=0
            for idx,prior in terms:
                v^=times(row[idx],previous[prior])
            output[mask]=v
        return output

    base=[0]*1024
    base[0]=1
    one=[extend(base,r,1) for r in range(10)]
    two={(i,j):extend(one[i],j,2) for i,j in combinations(range(10),2)}
    zeros=0
    count=0
    for c in combinations(range(10),3):
        start=extend(two[c[:2]],c[2],3)
        d_one=[extend(start,10+r,4) for r in range(8)]
        d_two={(i,j):extend(d_one[i],10+j,5) for i,j in combinations(range(9),2)}
        for a,b,t in combinations(range(10),3):
            final=extend(d_two[a,b],10+t,6)
            for mask,_ in indexed[6]:
                if final[mask]:
                    count+=1
                else:
                    zeros+=1
    return count,zeros

def main():
    outer=json.load(sys.stdin)["orientation"]
    c,d=chart_matrices(outer,WITNESSES[outer])
    checked,zeros=compute(tuple(zip(*c))+tuple(zip(*d)))
    assert checked==3024000 and zeros==0
    print(json.dumps({"orientation":outer,"nonzero_minors":checked,"zero_minors":zeros}))


if __name__=="__main__":
    main()
