"""The explicit four-group field code and its fixed tables."""
from functools import lru_cache
from itertools import combinations
from math import comb

M = 500
MASK = (1 << M) - 1
WIDTH = 50
CHUNKS = 10
SUBMASK = (1 << WIDTH) - 1
MOD_LOW = 1 | (1 << 1) | (1 << 2) | (1 << 9)
MODULUS = (1 << WIDTH) | MOD_LOW
THETA = 1 << WIDTH
LAMBDA = 7
EXTENSION_LOW = ((0,184),(1,108),(3,1))
QCOEFF = 1 << WIDTH
WITNESS = 0x1943d5008511e9019461d2f5cfb9bef6a69c36784b665ca1bf815d33d683f91b0ea077d4b578c5071547570e6904e2a692b76fdaed3624dbf196a3bea6267
THIRD = MASK // 3
CUBE_EXP = (2 * THIRD + 1) // 3
RADII = (99, 99, 101, 99, 101, 101)
SPLIT_RADIUS = 99
MAX_RADIUS = max(*RADII, SPLIT_RADIUS)
BYTES_PER_CHUNK = (WIDTH+7)//8
BYTES = CHUNKS * BYTES_PER_CHUNK
WORDS = CHUNKS
PAIRS = tuple(combinations(range(4), 2))
CHART_COSTS = (32000, 32000, 220000, 32000, 220000, 440000)
SPLIT_COST = 2200
CHOICES = tuple(combinations(range(CHUNKS), 5))
SPLIT_CHARTS = tuple((outer, known, rows) for outer in (0,1)
                     for known in CHOICES for rows in CHOICES)
ARRAY_LIST = (
    "Y0", "Y1", "Y2", "Y3", "KEY0", "KEY1", "V0", "V1", "V2", "V3",
    "BEST0", "BEST1", "X2", "X4", "T", "G", "H", "J", "CHECK", "C1", "C2",
    "INV", "INV2", "TEMP", "K", "TARGET", "ASROOT", "A0", "A1", "BASEA", "BASEB",
    "ROOT", "INVA", "KNOWN",
    "DELTA", "F1", "F2", "F3", "LIN", "QUAD", "ROW_TMP", "BASEINV",
) + tuple("MATRIX"+str(i) for i in range(CHUNKS)) + tuple(
    "COEF"+str(i) for i in range(CHUNKS)) + tuple("GAUSS"+str(i) for i in range(5))
BASE = {name: 48 + WORDS*i for i,name in enumerate(ARRAY_LIST)}
NEXT = 48 + WORDS * len(ARRAY_LIST)
LAST = NEXT + MAX_RADIUS + 1
INEXT = LAST + MAX_RADIUS + 1
ILAST = INEXT + MAX_RADIUS + 1
MULT = ILAST + MAX_RADIUS + 1
REDSHIFT1 = MULT + 16
REDSHIFT4 = REDSHIFT1 + 2
MAP_BASE = REDSHIFT4 + 16


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


def packed_words(a):
    yield from parts(a)


def cube(a):
    return product(a, square(a))


def frobenius(a, s):
    for _ in range(s):
        a = square(a)
    return a


def power(a, e):
    result = 1
    for bit in bin(e)[2:]:
        result = square(result)
        if bit == "1":
            result = product(result, a)
    return result


def geom(a, length, step):
    result = a
    k = 1
    for bit in bin(length)[3:]:
        result = product(result, frobenius(result, step * k))
        if bit == "1":
            result = product(result, frobenius(a, 2 * step * k))
        k = 2 * k + int(bit)
    return result


def inverse(a):
    assert a
    return square(geom(a, M - 1, 1))


def root_power(a):
    a7 = product(product(square(a), frobenius(a, 2)), a)
    t = frobenius(a7, 3)
    v = product(a, geom(t, 83, 6))
    return v


def cube_root(a):
    v = root_power(a)
    return v if cube(v) == a else None


@lru_cache(None)
def constant(tag):
    if tag == "lambda":
        return LAMBDA
    if tag == "lambda_inv":
        return inverse(LAMBDA)
    if tag == "lambda_minus_e":
        return power(LAMBDA, MASK - CUBE_EXP)
    if tag == "zeta":
        return power(LAMBDA, THIRD)
    if tag == "zeta2":
        return square(constant("zeta"))
    if isinstance(tag,str) and tag.startswith("theta"):
        return power(THETA, int(tag[5:]))
    raise ValueError("unknown field constant")


@lru_cache(None)
def AS_pivots():
    pivots = {}
    for i in range(M):
        preimage = 1 << i
        image = square(preimage) ^ preimage
        while image:
            pivot = image.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (image, preimage)
                break
            r, x = pivots[pivot]
            image ^= r
            preimage ^= x
    return tuple((i, *pivots[i]) for i in sorted(pivots, reverse=True))


def AS_map(y):
    x = 0
    for pivot, image, preimage in AS_pivots():
        if (y >> pivot) & 1:
            y ^= image
            x ^= preimage
    return x


def cube_roots(x):
    base = cube_root(x)
    if base is None:
        return ()
    return (base, product(constant("zeta"), base), product(constant("zeta2"), base))


def codeword(a, b):
    return a, b, product(a, b), cube(a) ^ product(LAMBDA, cube(b))


def encode(message):
    coordinates = codeword(message & MASK, message >> M)
    return sum(value << (M * i) for i, value in enumerate(coordinates))


def chart(pair, x, y):
    """Chart candidates passing the good-message test."""
    if pair == (2, 3) and not y:
        return ()
    if not x:
        return ()
    return tuple((a,b) for a,b in raw_chart(pair,x,y) if a and b)


def raw_chart(pair, x, y):
    """Host semantic counterpart of the bounded field chart macros."""
    if pair == (0, 1):
        return ((x, y),) if x and y else ()
    if not x:
        return ()
    if pair == (0, 2):
        return ((x, product(y, inverse(x))),) if y else ()
    if pair == (1, 2):
        return ((product(y, inverse(x)), x),) if y else ()
    if pair == (0, 3):
        target = product(y ^ cube(x), constant("lambda_inv"))
        return tuple((x, b) for b in cube_roots(target)) if target else ()
    if pair == (1, 3):
        target = y ^ product(LAMBDA, cube(x))
        return tuple((a, x) for a in cube_roots(target)) if target else ()
    assert pair == (2, 3)
    c, d = x, y
    if not d:
        return ()
    k = product(LAMBDA, cube(c))
    z = product(k, square(inverse(d)))
    w = AS_map(z)
    if square(w) ^ w != z:
        return ()
    aa = product(d, w)
    root = root_power(aa)
    test = cube(root)
    if test == aa:
        a, b = root, product(c, inverse(root))
    elif test == product(constant("zeta2"), aa):
        b = product(root, constant("lambda_minus_e"))
        a = product(c, inverse(b))
    else:
        return ()
    return ((a, b),
            (product(constant("zeta"), a), product(constant("zeta2"), b)),
            (product(constant("zeta2"), a), product(constant("zeta"), b)))


def map_tags():
    shifts = {1, 2, 3}
    for length, step in ((83, 6), (M - 1, 1)):
        k = 1
        for bit in bin(length)[3:]:
            shifts.add(step * k)
            if bit == "1":
                shifts.add(2 * step * k)
            k = 2 * k + int(bit)
    return tuple(sorted(shifts)) + ("lambda", "lambda_inv", "lambda_minus_e", "zeta", "zeta2", "AS") + tuple("theta"+str(i) for i in sorted(set(range(1,10))|set(range(2,19,2))) )


def map_columns(tag):
    for chunk in range(CHUNKS):
        for i in range(8*(BYTES//CHUNKS)):
            if i >= WIDTH:
                yield 0
                continue
            e = 1 << (chunk*WIDTH+i)
            if tag == "AS":
                yield AS_map(e)
            elif isinstance(tag, int):
                yield frobenius(e, tag)
            else:
                yield product(constant(tag),e)


def map_words(tag):
    cols=iter(map_columns(tag))
    for _ in range(BYTES):
        vals=[0]
        for _ in range(8):
            v=next(cols)
            vals.extend(x^v for x in vals[:])
        for value in vals:
            yield from packed_words(value)


def multiplication_columns(v):
    return tuple(parts(product(v,1 << (WIDTH*i))) for i in range(CHUNKS))


def complement(indices):
    return tuple(i for i in range(CHUNKS) if i not in indices)


def coefficient_matrix(v,known,rows):
    unknown=complement(known)
    cols=multiplication_columns(v)
    return tuple(tuple(cols[j][i] for j in unknown) for i in rows)


def inverse_matrix(matrix):
    n=len(matrix)
    augmented=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if augmented[i][j]),None)
        if pivot is None:
            return None
        augmented[j],augmented[pivot]=augmented[pivot],augmented[j]
        inv=base_power(augmented[j][j],SUBMASK-1)
        augmented[j]=[base_product(inv,x) for x in augmented[j]]
        for i in range(n):
            if i==j:
                continue
            scale=augmented[i][j]
            augmented[i]=[x^base_product(scale,y) for x,y in zip(augmented[i],augmented[j])]
    return tuple(tuple(row[n:]) for row in augmented)


def split_map(v,known,rows):
    M=coefficient_matrix(v,known,rows)
    inv=inverse_matrix(M)
    if inv is None:
        return None
    cols=multiplication_columns(v)
    out=[]
    for row in inv:
        products=[]
        for j in known:
            value=0
            for coefficient,i in zip(row,rows):
                value^=base_product(coefficient,cols[j][i])
            products.append(value)
        out.append(tuple(products)+row)
    return tuple(out)


def descriptor(outer,known,rows):
    return (outer,*known,*complement(known),*rows)

def split_recover(v, known, rows, key):
    coeff=split_map(v,known,rows)
    if coeff is None:
        return None
    u=[0]*CHUNKS
    for at,val in zip(known,key):
        u[at]=val
    for i,at in enumerate(complement(known)):
        for c,val in zip(coeff[i],key):
            u[at]^=base_product(c,val)
    uu=join(u)
    return uu,product(v,uu)


def directions(v,known,rows):
    coeff=split_map(v,known,rows)
    assert coeff is not None
    for lane in range(CHUNKS):
        for i in range(WIDTH):
            e=1 << i
            u=[0]*CHUNKS
            if lane < len(known):
                u[known[lane]]=e
            for pos,at in enumerate(complement(known)):
                u[at]=base_product(coeff[pos][lane],e)
            uu=join(u)
            yield uu,product(v,uu)


def derivative_columns(alpha,e):
    lin=product(alpha,square(e))
    quad=product(alpha,e)
    for lane in range(CHUNKS):
        l=lin if lane == 0 else product(constant(f"theta{lane}"),lin)
        q=quad if lane == 0 else product(constant(f"theta{2*lane}"),quad)
        for _ in range(WIDTH):
            yield l^q
            l=product(l,2)
            q=product(q,4)



AFFINE_WIDTHS = (13,13,13,11)
EVAL_WIDTHS = (20,20,10)
EVAL_POINTS = tuple(range(19))


def affine_entry(alpha,e,slice_index,value):
    lane,part=divmod(slice_index,len(AFFINE_WIDTHS))
    begin=lane*WIDTH+sum(AFFINE_WIDTHS[:part])
    cols=tuple(derivative_columns(alpha,e))
    out=product(alpha,cube(e)) if slice_index==0 else 0
    for bit in range(AFFINE_WIDTHS[part]):
        if (value >> bit)&1:
            out^=cols[begin+bit]
    return out


@lru_cache(None)
def evaluation_columns():
    return tuple(tuple(base_power(x,j) if j else 1 for x in EVAL_POINTS)
                 for j in range(CHUNKS))


@lru_cache(None)
def interpolation_columns():
    result=[]
    for x in EVAL_POINTS:
        terms=[1]
        scale=1
        for y in EVAL_POINTS:
            if y==x:
                continue
            updated=[0]*(len(terms)+1)
            for j,a in enumerate(terms):
                updated[j]^=base_product(a,y)
                updated[j+1]^=a
            terms=updated
            scale=base_product(scale,x^y)
        inverse=base_power(scale,SUBMASK-1)
        result.append(parts(reduce_extension([base_product(t,inverse) for t in terms])))
    return tuple(result)


def linear_slice_entry(columns,position,width,lane,key):
    assert 0 <= key < 1 << width
    value=key << position
    return tuple(base_product(value,c) for c in columns[lane])


def evaluate_product(a,b):
    columns=evaluation_columns()
    evals=[[0]*len(EVAL_POINTS) for _ in range(2)]
    for values,field in zip(evals,(a,b)):
        for coefficient,row in zip(parts(field),columns):
            for j,c in enumerate(row):
                values[j]^=base_product(coefficient,c)
    product_at_points=tuple(base_product(x,y) for x,y in zip(*evals))
    coefficients=[0]*CHUNKS
    for value,row in zip(product_at_points,interpolation_columns()):
        for j,c in enumerate(row):
            coefficients[j]^=base_product(value,c)
    return join(coefficients)


def nonzero_minor_count():
    matrix=tuple(zip(*multiplication_columns(WITNESS)))
    maps={}
    for v in set(x for row in matrix for x in row):
        columns=[v]
        for _ in range(WIDTH-1):
            columns.append(remainder(columns[-1]<<1))
        tables=[]
        for j in range(BYTES_PER_CHUNK):
            values=[0]
            for e in range(8):
                i=8*j+e
                x=columns[i] if i<WIDTH else 0
                values.extend(y^x for y in values[:])
            tables.append(values)
        maps[v]=tables

    def multiply(tables,x):
        value=0
        for table in tables:
            value^=table[x&255]
            x>>=8
        return value

    previous={((),()):1}
    for k in range(1,6):
        sets=tuple(combinations(range(CHUNKS),k))
        result={}
        for cols in sets:
            omitted=[cols[:i]+cols[i+1:] for i in range(k)]
            for rows in sets:
                literal=matrix[rows[-1]]
                value=0
                for col,prior_cols in zip(cols,omitted):
                    value^=multiply(maps[literal[col]],previous[rows[:-1],prior_cols])
                result[rows,cols]=value
        previous=result
    return sum(value!=0 for value in previous.values())


def reduction_table(s):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << s))


CHART_META = MAP_BASE + len(map_tags()) * BYTES * 256 * WORDS
DIRECTIONS = CHART_META + len(SPLIT_CHARTS)*16
SAVED = DIRECTIONS + M * 2 * WORDS
AFFINE = SAVED + (MAX_RADIUS+1) * 3 * WORDS
AFFINE_STRIDE = CHUNKS * sum(1<<b for b in AFFINE_WIDTHS) * WORDS
EVAL_A = AFFINE + M*AFFINE_STRIDE
EVAL_B = EVAL_A + len(EVAL_POINTS)
EVAL_PRODUCT = EVAL_B + len(EVAL_POINTS)
EVAL_TABLE = EVAL_PRODUCT + len(EVAL_POINTS)
LINEAR_ROWS = CHUNKS * sum(1<<b for b in EVAL_WIDTHS) * len(EVAL_POINTS)
INTERP_TABLE = EVAL_TABLE + LINEAR_ROWS
MEMORY_END = INTERP_TABLE + LINEAR_ROWS


def W_bound():
    sizes=tuple(sum(comb(2*M,i) for i in range(r+1)) for r in RADII)
    outer=tuple(sum(comb(M,i) for i in range(r+1)) for r in RADII)
    ns=sum(comb(2*M,i) for i in range(SPLIT_RADIUS+1))
    os=sum(comb(M,i) for i in range(SPLIT_RADIUS+1))
    ordinary=sum((c+40)*n+20000000062*o for c,n,o in zip(CHART_COSTS,sizes,outer))
    split=len(SPLIT_CHARTS)*((SPLIT_COST+40)*ns+20000000062*os+200)
    return sizes,outer,ns,os,ordinary+split+2000
