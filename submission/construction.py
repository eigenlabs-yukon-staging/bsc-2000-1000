"""Complete fixed systematic block-code construction; no external data."""
from array import array

P48 = [10062753, 14743013, 14914939, 9919255, 9168521, 4546510, 7708315, 3300079, 5800557, 12747841, 2443661, 9067217, 3382048, 3666711, 13970140, 107303, 12645428, 16635897, 11525756, 823522, 13926933, 5833403, 9817670, 7719323]
P32 = [32767,33383,33963,35123,37069,41301,49561,38414,43542,52250,61468,38624,43872,52640,61888,65024]
T48=48
T32=T48+(1<<24)
LUT48=T32+(1<<16)
LUT32=LUT48+6*256
END=LUT32+4*256

def encode_block(data, parity):
    p=0
    for i,col in enumerate(parity):
        if (data>>i)&1:p^=col
    return data | (p<<len(parity))

def encode(message):
    out=0
    for b in range(41):
        out |= encode_block((message>>(24*b))&((1<<24)-1),P48)<<(48*b)
    out |= encode_block((message>>984)&65535,P32)<<1968
    return out

def syndrome_luts(parity):
    k=len(parity); cols=parity+[1<<i for i in range(k)]
    tables=array('I')
    for chunk in range(2*k//8):
        for value in range(256):
            s=0
            for bit in range(8):
                if value>>bit&1:s^=cols[8*chunk+bit]
            tables.append(s)
    return tables

def leaders(parity):
    """BFS from syndrome zero; neighbours in bit order; first discovery wins.

    Returns the low-k-bit information correction of each chosen full leader,
    and the exact histogram of its minimum weights. Every table entry is built.
    """
    k=len(parity); size=1<<k
    columns=parity+[1<<i for i in range(k)]
    toggles=[1<<i for i in range(k)]+[0]*k
    edges=list(zip(columns,toggles))
    dist=bytearray([255])*size;dist[0]=0
    correction=array('I',[0])*size
    queue=array('I',[0]);append=queue.append
    pos=0
    while len(queue)<size:
        s=queue[pos];pos+=1
        d=dist[s]+1;v0=correction[s]
        for col,toggle in edges:
            v=s^col
            if dist[v]==255:
                dist[v]=d;correction[v]=v0^toggle;append(v)
    hist=[dist.count(i) for i in range(max(dist)+1)]
    assert sum(hist)==size and len(queue)==size
    return correction,hist
