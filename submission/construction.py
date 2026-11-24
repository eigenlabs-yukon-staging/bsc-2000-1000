"""Fixed field arithmetic and six-million-chart bilinear decoder."""
from functools import lru_cache
from itertools import combinations
from math import comb

WIDTH=50
CHUNKS=10
M=WIDTH*CHUNKS
MASK=(1<<M)-1
SUBMASK=(1<<WIDTH)-1
MOD_LOW=1|(1<<1)|(1<<2)|(1<<9)
MODULUS=(1<<WIDTH)|MOD_LOW
THETA=1<<WIDTH
EXTENSION_LOW=((0,184),(1,108),(3,1))
RADIUS=86
BYTES=70
KNOWN_SETS=tuple(combinations(range(10),4))
PARITY_SETS=tuple(combinations(range(10),3))
ORIENTATION_CHARTS=len(KNOWN_SETS)*len(PARITY_SETS)**2
SPLIT_COUNT=2*ORIENTATION_CHARTS

ARRAYS=(
    "Y0","Y1","Y2","Y3","KEY0","KEY1","V0","V1","V2","V3","BEST0","BEST1",
    "VROOT","TEMP","PROD0","PROD1","PRODTEMP","RED_SLOT","DELTA","DC","DD",
    "ROW_TMP","BASEINV","KEY_TRANSFORM",
)+tuple("MC"+str(i) for i in range(10))+tuple(
    "MD"+str(i) for i in range(10))+tuple(
    "COEF"+str(i) for i in range(10))+tuple("GAUSS"+str(i) for i in range(8))
BASE={name:48+10*i for i,name in enumerate(ARRAYS)}
NEXT=48+10*len(ARRAYS)
LAST=NEXT+RADIUS+1
INEXT=LAST+RADIUS+1
ILAST=INEXT+RADIUS+1
MULT=ILAST+RADIUS+1
REDSHIFT1=MULT+16
REDSHIFT4=REDSHIFT1+2
GAMMA_MAP=REDSHIFT4+16
MAP_BASE=GAMMA_MAP+2*7*256
MAP_TAGS=("sqrt",)+tuple(f"theta{i}" for i in range(1,10))+tuple(
    f"root_theta{i}" for i in range(1,10))
CHART_META=MAP_BASE+len(MAP_TAGS)*BYTES*256*10
DIRECTIONS=CHART_META+SPLIT_COUNT*17
MEMORY_END=DIRECTIONS+M*30


def remainder(a, b=MODULUS):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def carryless(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def base_product(a, b):
    return remainder(carryless(a, b))


def base_power(a, e):
    out = 1
    for bit in bin(e)[2:]:
        out = base_product(out, out)
        if bit == "1":
            out = base_product(out, a)
    return out


def parts(a):
    return tuple((a >> (WIDTH*i)) & SUBMASK for i in range(CHUNKS))


def join(values):
    return sum(v << (WIDTH*i) for i,v in enumerate(values))


def reduce_extension(terms):
    for k in range(2*CHUNKS-2,CHUNKS-1,-1):
        if not terms[k]:
            continue
        for j,c in EXTENSION_LOW:
            terms[k-CHUNKS+j]^=base_product(terms[k],c)
    return join(terms[:CHUNKS])


def product(a,b):
    terms=[0]*(2*CHUNKS-1)
    for i,x in enumerate(parts(a)):
        if x:
            for j,y in enumerate(parts(b)):
                if y:
                    terms[i+j]^=base_product(x,y)
    return reduce_extension(terms)


def square(a):
    terms=[0]*(2*CHUNKS-1)
    for i,x in enumerate(parts(a)):
        terms[2*i]=base_product(x,x)
    return reduce_extension(terms)


def field_words(a):
    yield from parts(a)



@lru_cache(maxsize=128)
def sqrt(a):
    for _ in range(M-1):
        a=square(a)
    return a


def codeword(a,b):
    return a,b,product(a,b),product(sqrt(a),b)


def encode(message):
    return sum(value<<(M*i) for i,value in enumerate(codeword(message&MASK,message>>M)))


@lru_cache(None)
def scalar(tag):
    if tag.startswith("root_"):
        return sqrt(scalar(tag[5:]))
    return 1<<(WIDTH*int(tag[5:]))


def table_columns(tag):
    for i in range(CHUNKS):
        for j in range(56):
            if j>=WIDTH:
                yield 0
            else:
                value=1<<(WIDTH*i+j)
                yield sqrt(value) if tag=="sqrt" else product(value,scalar(tag))


def table_words(tag):
    cols=iter(table_columns(tag))
    for _ in range(BYTES):
        row=[0]
        for _ in range(8):
            x=next(cols)
            row.extend(y^x for y in row[:])
        for x in row:
            yield from field_words(x)


def gamma_entry(index,j,value):
    constant=EXTENSION_LOW[index][1]
    return base_product(constant,remainder(value<<(8*j)))


def reduction_table(s):
    return tuple(carryless(x,MOD_LOW) for x in range(1<<s))


def chart_at(ordinal):
    outer,relative=divmod(ordinal,ORIENTATION_CHARTS)
    key,p=divmod(relative,len(PARITY_SETS)**2)
    c,d=divmod(p,len(PARITY_SETS))
    return outer,KNOWN_SETS[key],PARITY_SETS[c],PARITY_SETS[d]


def unknown_indices(known):
    return tuple(i for i in range(CHUNKS) if i not in known)


def descriptor(ordinal):
    outer,known,cr,dr=chart_at(ordinal)
    return (outer,*known,*unknown_indices(known),*cr,*dr)


@lru_cache(maxsize=8)
def chart_matrices(outer,v):
    c=tuple(parts(product(v,1<<(WIDTH*j))) for j in range(CHUNKS))
    if outer==0:
        root=sqrt(v)
        d=tuple(parts(product(root,1<<(WIDTH*j))) for j in range(CHUNKS))
    else:
        d=tuple(tuple(base_product(x,x)
                      for x in parts(product(v,1 if j==0 else scalar(f"root_theta{j}"))))
                for j in range(CHUNKS))
    return c,d


def inverse_matrix(m):
    n=len(m)
    rows=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(m)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if rows[i][j]),None)
        if pivot is None:
            return None
        rows[j],rows[pivot]=rows[pivot],rows[j]
        factor=base_power(rows[j][j],SUBMASK-1)
        rows[j]=[base_product(factor,x) for x in rows[j]]
        for i in range(n):
            if i!=j:
                scale=rows[i][j]
                rows[i]=[x^base_product(scale,y) for x,y in zip(rows[i],rows[j])]
    return tuple(tuple(row[n:]) for row in rows)


def dot(a,b):
    result=0
    for x,y in zip(a,b):
        result^=base_product(x,y)
    return result


@lru_cache(maxsize=32)
def recovery(chart,v):
    outer,known,cr,dr=chart
    unknown=unknown_indices(known)
    c,d=chart_matrices(outer,v)
    rows=tuple(tuple(cols[i][j] for i in range(CHUNKS))
               for cols,inds in ((c,cr),(d,dr)) for j in inds)
    unknown_matrix=tuple(tuple(row[i] for i in unknown) for row in rows)
    inverse=inverse_matrix(unknown_matrix)
    if inverse is None:
        return None
    known_columns=tuple(tuple(row[i] for row in rows) for i in known)
    return unknown,tuple(tuple(dot(inv,col) for col in known_columns)+inv
                         for inv in inverse)


def recover(chart,v,key):
    outer,known,_,_=chart
    result=recovery(chart,v)
    if result is None:
        return None
    unknown,coeff=result
    if outer:
        key=key[:7]+tuple(base_product(x,x) for x in key[7:])
    out=[0]*CHUNKS
    for at,val in zip(known,key):
        out[at]=val
    for at,col in zip(unknown,coeff):
        out[at]=dot(col,key)
    u=join(out)
    return codeword(v,u) if outer==0 else codeword(u,v)


def direction(chart,v,system,i):
    outer,known,_,_=chart
    unknown,coeff=system
    lane,bit=divmod(i,WIDTH)
    f=1<<bit
    transformed=base_product(f,f) if outer and lane>=7 else f
    row=[0]*CHUNKS
    if lane<4:
        row[known[lane]]=f
    for at,co in zip(unknown,coeff):
        row[at]=base_product(co[lane],transformed)
    e=join(row)
    dc=product(v,e)
    dd=product(sqrt(v),e) if outer==0 else product(v,sqrt(e))
    return e,dc,dd


def operation_bound():
    n=sum(comb(1000,i) for i in range(RADIUS+1))
    o=sum(comb(500,i) for i in range(RADIUS+1))
    return n,o,SPLIT_COUNT*(790*n+50000062*o+200)+100000
