"""Reproduce the fixed local spectra and the exact polynomial-tag certificate."""
import base64
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import sys
import zlib

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from collections import Counter
from functools import lru_cache
from construction import (
    PRIVATE_BITS, PARITY_BITS, SIDE_BITS, TAG_BITS, SLOTS, COLUMNS,
    SIDES, DIMENSION, COMPONENT_BITS, PAYLOAD_BITS, COMPONENTS, RADIUS,
    MODULUS, EVALUATION_POINTS, ANCHORS, COEFF_BASE, COEFF_END, DESCRIPTOR_BASE,
    DESCRIPTOR_END, ROW_STRIDE, PARITY_BASE, TAG_BASE, LOW_SPLIT, HIGH_SPLIT,
    LIST_BASE, TAG_SCRATCH, TAG_SCRATCH_END, LEADER_BASE, gf_mul, gf_pow,
    coefficients, coefficient_key, table_word, descriptor_base, list_address,
    displacement, split_basis, split_syndrome, syndrome_columns, linear_map, schedule,
)

EXACT_DEPTH = 11
RANK_CAP = 64
SECTOR_BOUNDARIES = (0, 1100, 2200, 3400, 4096)


def gray_order(columns):
    values = [0]
    for column in columns:
        values.extend(value ^ column for value in reversed(values))
    return values


def read_histograms():
    streams = []
    for i in range(4):
        data = json.loads((HERE / f"spectra{i}.json").read_text())
        raw_block = zlib.decompress(base64.b64decode(data["packed_spectra"], validate=True))
        assert len(raw_block) == (SLOTS // 4) * (SIDE_BITS + 1) * 3
        streams.append(raw_block)
    raw = b"".join(streams)
    rows = []
    for tag in range(SLOTS):
        start = (SIDE_BITS + 1) * 3 * tag
        rows.append([int.from_bytes(raw[start+3*i:start+3*i+3], "little") for i in range(SIDE_BITS+1)])
        assert sum(rows[-1]) == 1 << PRIVATE_BITS
    return rows

def check_sector(sector, hists):
    # A word lists the private bits first, followed by the parity bits.
    basis = tuple((1 << i) | (col << PRIVATE_BITS) for i, col in enumerate(COLUMNS[:PRIVATE_BITS]))
    words = gray_order(basis)
    for tag in range(SECTOR_BOUNDARIES[sector], SECTOR_BOUNDARIES[sector + 1]):
        shift = displacement(tag) << PRIVATE_BITS
        hist = [0] * (SIDE_BITS + 1)
        for word in words:
            hist[(word ^ shift).bit_count()] += 1
        assert hist == hists[tag]


def coverage(volumes, pair_limits):
    cumulative = []
    previous = 0
    for volume, pairs in zip(volumes, pair_limits, strict=True):
        k = 1 + pairs // volume
        denom = k * (k + 1)
        lower = (2 * k * volume - pairs + denom - 1) // denom
        cauchy = (volume * volume + volume + pairs - 1) // (volume + pairs)
        previous = max(previous, lower, cauchy)
        cumulative.append(previous)
    return cumulative


RANK_QUOTIENT_CHECKS = (8193,65553)
RANK_QUOTIENT_PIVOTS = (0,4)


def quotient_sector_coordinate(q):
    inner=sum(((q>>bit)&1)<<i for i,bit in enumerate(
        b for b in range(PRIVATE_BITS) if b not in RANK_QUOTIENT_PIVOTS))
    sector=sum(((q&row).bit_count()&1)<<i for i,row in enumerate(RANK_QUOTIENT_CHECKS))
    return inner | (sector<<(PRIVATE_BITS-2))


def exact_rank_sector(sector):
    from functools import reduce
    from itertools import combinations
    from operator import xor

    tag_mask=(1<<TAG_BITS)-1
    row_bits=PRIVATE_BITS-2
    size=1<<row_bits
    local_mask=(1<<(row_bits+TAG_BITS))-1
    columns=[]
    for word in syndrome_columns():
        q=split_syndrome(word)
        columns.append((quotient_sector_coordinate(q>>TAG_BITS)<<TAG_BITS) | (q&tag_mask))
    p,q=RANK_QUOTIENT_PIVOTS
    a,b=RANK_QUOTIENT_CHECKS
    assert ((a>>p)&1)*((b>>q)&1) ^ ((a>>q)&1)*((b>>p)&1) == 1
    # Normalize tag coordinates within each quotient using an
    # F2-linear function of the quotient.
    pivots={}
    for c in columns:
        if c>local_mask:
            continue
        q,tag=c>>TAG_BITS,c&tag_mask
        while q:
            bit=q.bit_length()-1
            if bit not in pivots:
                pivots[bit]=q,tag
                break
            old,shift=pivots[bit]
            q,tag=q^old,tag^shift
    assert len(pivots)==row_bits

    def normalize(c):
        q=(c>>TAG_BITS)&(size-1)
        tag=0
        while q:
            old,shift=pivots[q.bit_length()-1]
            q,tag=q^old,tag^shift
        return c^tag

    columns=list(map(normalize,columns))
    assert sum(c<=local_mask and (c&tag_mask)==0 for c in columns)>=row_bits
    single=tuple(1<<i for i in range(1<<TAG_BITS))
    full=(1<<(1<<TAG_BITS))-1
    shifts=[]
    for bit in range(TAG_BITS):
        step=1<<bit
        lo=full//((1<<step)+1)
        shifts.append((step,lo,lo<<step))
    preserving=tuple(sorted((c for c in columns if c<=local_mask),
                            key=lambda c:(c&tag_mask).bit_count()))
    other=tuple(c for c in columns if c>local_mask)
    assert len(preserving)==36 and len(other)==14
    translations=[
        (c>>TAG_BITS,tuple(shifts[i] for i in range(TAG_BITS) if (c>>i)&1))
        for c in preserving
    ]
    rows=[0]*size
    previous=[0]*size
    history=[]
    base=sector<<(row_bits+TAG_BITS)
    top=base+local_mask
    for w in range(EXACT_DEPTH+1):
        if w<=5:
            direct=columns
        else:
            use=min(len(preserving),26+3*(w-6))
            next_rows=rows.copy()
            if sum(previous)*8>=7*size*(1<<TAG_BITS):
                for j,union in enumerate(rows):
                    if union==full:
                        continue
                    for high,swaps in translations[:use]:
                        row=rows[j^high]
                        for step,lo,hi in swaps:
                            row=((row&lo)<<step)|((row&hi)>>step)
                        union|=row
                        if union==full:
                            break
                    next_rows[j]=union
            else:
                for high,swaps in translations[:use]:
                    if not swaps:
                        for j,row in enumerate(rows):
                            next_rows[j^high]|=row
                    else:
                        for j,row in enumerate(rows):
                            if row:
                                for step,lo,hi in swaps:
                                    row=((row&lo)<<step)|((row&hi)>>step)
                                next_rows[j^high]|=row
            rows=next_rows
            direct=preserving[use:]+other
        for group in combinations(direct,w):
            s=reduce(xor,group,0)
            if base<=s<=top:
                rows[(s&local_mask)>>TAG_BITS]|=single[s&tag_mask]
        delta=[0]*(RANK_CAP+1)
        above=0
        for j,row in enumerate(rows):
            now=row.bit_count()
            old=previous[j]
            assert old<=now
            first,last=min(old,RANK_CAP),min(now,RANK_CAP)
            delta[first]+=1
            delta[last]-=1
            above+=now-old-(last-first)
            previous[j]=now
        ranks=[]
        running=0
        for x in delta[:-1]:
            running+=x
            ranks.append(running)
        ranks.append(above)
        history.append(ranks)
    return history


def rank_witness():
    data=json.loads((HERE/'ranks.json').read_text())
    assert len(data)==4 and RADIUS<=RANK_CAP
    assert all(len(sector)==EXACT_DEPTH+1 for sector in data)
    assert all(len(row)==RANK_CAP+1 and all(type(x) is int and x>=0 for x in row)
               for sector in data for row in sector)
    assert all(sum(map(sum,sector))==1<<(PARITY_BITS-2) for sector in data)
    total=[[sum(counts) for counts in zip(*rows,strict=True)]
           for rows in zip(*data,strict=True)]
    return data,total


def check_rank_sector(sector, witness):
    assert exact_rank_sector(sector)==witness


PROFILE_SECTORS=((5,),(6,),(7,),(8,))
EXPECTED_W=611362

GLOBAL_CUTOFF=54
SCALE=1<<24
STRIDE=512
DIGIT_MASK=(1<<STRIDE)-1
POLY_MODULUS=1<<(STRIDE*(2*GLOBAL_CUTOFF+1))


def ceildiv(a,b):
    return -(-a//b)


@lru_cache(None)
def root_ceil(value,exponent):
    lo=0
    hi=1<<ceildiv(value.bit_length(),exponent)
    while lo<hi:
        mid=(lo+hi)//2
        if mid**exponent>=value:
            hi=mid
        else:
            lo=mid+1
    return lo


def coefficient(poly,weight):
    return (poly>>(STRIDE*weight))&DIGIT_MASK


def checked_pack(counts):
    assert all(0<=x<=DIGIT_MASK for x in counts)
    return sum(x<<(STRIDE*i) for i,x in enumerate(counts))


def exact_five(kernels,history):
    from itertools import combinations

    roots=EVALUATION_POINTS[:3]
    field=set(EVALUATION_POINTS)
    assert len(field)==8 and 0 in field and 1 in field
    assert all(x^y in field and gf_mul(x,y) in field for x in field for y in field)
    orbit={tuple(sorted(gf_mul(a,x)^b for x in roots))
           for a in field if a for b in field}
    assert orbit==set(combinations(EVALUATION_POINTS,3))
    values=[]
    for point in EVALUATION_POINTS[3:]:
        result=1
        for root in roots:
            result=gf_mul(result,point^root)
        assert result
        values.append(result)

    stride=288
    mask=(1<<stride)-1
    modulus=1<<(stride*(2*GLOBAL_CUTOFF+1))
    totals=[sum(row) for row in history]
    assert all(sum(row)<=1<<PARITY_BITS for row in kernels)
    assert (SLOTS-1)*sum(totals)**SIDES<1<<stride
    encoded=[0]+[sum(x<<(stride*w) for w,x in enumerate(row)) for row in kernels]
    outer=0
    for scalar in range(1,SLOTS):
        product=1
        for value in values:
            product=(product*encoded[gf_mul(scalar,value)])%modulus
        outer+=product
    common=sum(total<<(stride*2*w) for w,total in enumerate(totals))
    packed=(outer*pow(common,3,modulus))%modulus
    result=[]
    running=0
    for w in range(2*GLOBAL_CUTOFF+1):
        running+=(packed>>(stride*w))&mask
        if w%2==0:
            result.append(running)
    return result


def local_pair_shells(hists, history):
    totals=[sum(row) for row in history]
    assert sum(totals)==1<<PARITY_BITS
    largest=len(totals)-1
    accumulated=[0]
    for value in totals:
        accumulated.append(accumulated[-1]+value)

    def marginal(w):
        return accumulated[min(max(w+1,0),largest+1)]

    thresholds=[min(marginal(k)+marginal(w-k-1)
                    for k in range(-1,largest+1)) for w in range(2*largest+1)]
    factors=[]
    for u in range(largest+1):
        for v in range(u,largest+1):
            pairs=[]
            for d in range(max(1,v-u),u+v+1):
                a=(d+u-v)//2
                b=(u+v-d)//2
                if (d+u-v)%2==0 and 0<=a<=d and 0<=b<=SIDE_BITS-d:
                    pairs.append((d,comb(d,a)*comb(SIDE_BITS-d,b)))
            factors.append((u,v,tuple(pairs)))

    result=[]
    for counts in hists[1:]:
        bounds=[[0]*(largest+1) for _ in totals]
        for u,v,pairs in factors:
            ambient=sum(counts[d]*value for d,value in pairs)
            bounds[u][v]=min(ambient,totals[v])
            bounds[v][u]=min(ambient,totals[u])

        for row in bounds:
            for v in range(1,largest+1):
                row[v]+=row[v-1]
        previous=0
        shells=[]
        for w,threshold in enumerate(thresholds):
            by_first=sum(min(totals[u],bounds[u][min(w-u,largest)])
                         for u in range(min(w,largest)+1))
            cap=min(threshold,by_first)
            assert cap>=previous
            shells.append(cap-previous)
            previous=cap
        result.append(shells)
    return result


def profile(changes, kernels, history):
    if changes==5:
        return exact_five(kernels,history)
    freedom=changes-DIMENSION
    norms=[]
    for w in range(2*len(history)-1):
        total=sum(root_ceil(SCALE**freedom*row[w]**changes,freedom)
                  for row in kernels if row[w])
        norms.append(root_ceil(SCALE**(changes-freedom)*total**freedom,changes))

    common=[0]*(2*len(history)-1)
    for w,row in enumerate(history):
        common[2*w]=sum(row)
    assert sum(norms)**changes*sum(common)**(SIDES-changes)<=DIGIT_MASK
    packed=(pow(checked_pack(norms),changes,POLY_MODULUS)*
            pow(checked_pack(common),SIDES-changes,POLY_MODULUS))%POLY_MODULUS
    running=0
    result=[]
    for d in range(2*GLOBAL_CUTOFF+1):
        running+=coefficient(packed,d)
        if d%2==0:
            result.append(ceildiv(running,SCALE**changes))
    return result


def profile_witness():
    data=json.loads((HERE/'profiles.json').read_text())
    assert len(data)==SIDES-DIMENSION
    assert all(len(row)==GLOBAL_CUTOFF+1 for row in data)
    assert all(isinstance(x,str) and x.isdecimal() for row in data for x in row)
    return {changes:list(map(int,row))
            for changes,row in enumerate(data,DIMENSION+1)}


def component_numerator(profiles,history):
    totals=[sum(row) for row in history]
    assert sum(totals)==1<<PARITY_BITS
    admitted=[sum(row[:RADIUS]) for row in history]
    rejected=[total-a for total,a in zip(totals,admitted,strict=True)]
    assert sum(admitted)==RADIUS*(1<<PRIVATE_BITS)
    pairs=[sum(comb(SIDES,c)*profiles[c][t] for c in profiles)
           for t in range(GLOBAL_CUTOFF+1)]
    a,u=checked_pack(admitted),checked_pack(rejected)
    assert sum(totals)**SIDES<=DIGIT_MASK
    packed=sum(comb(SIDES,i)*pow(a,i,POLY_MODULUS)*pow(u,SIDES-i,POLY_MODULUS)
               for i in range(DIMENSION,SIDES+1))%POLY_MODULUS
    volumes=[]
    cumulative=0
    for t in range(GLOBAL_CUTOFF+1):
        cumulative+=coefficient(packed,t)
        volumes.append(cumulative)
    occupied=coverage(volumes,pairs)
    assert occupied[-1]<=1<<(COMPONENT_BITS-PAYLOAD_BITS)
    numerator=sum((x-(occupied[t-1] if t else 0))*9**(COMPONENT_BITS-t)
                  for t,x in enumerate(occupied))
    return numerator,occupied


def polynomial_gcd(a,b):
    while b:
        while a.bit_length()>=b.bit_length():
            a^=b<<(a.bit_length()-b.bit_length())
        a,b=b,a
    return a


def check_program(exponent):
    assert MODULUS.bit_length()==TAG_BITS+1
    assert gf_pow(2,SLOTS)==2
    for p in (2,3):
        assert polynomial_gcd(MODULUS,gf_pow(2,1<<(TAG_BITS//p))^2)==1
    assert len(set(EVALUATION_POINTS))==SIDES
    assert COMPONENT_BITS*COMPONENTS==2000 and PAYLOAD_BITS*COMPONENTS==1000
    assert len(ANCHORS)==comb(SIDES,DIMENSION)==70
    for g,anchors in enumerate(ANCHORS):
        matrix=coefficients(anchors)
        for i in range(DIMENSION):
            assert tuple(matrix[i][s] for s in anchors)==tuple(int(i==j) for j in range(DIMENSION))
            for value in (0,1,2,4095):
                other=[s for s in range(SIDES) if s not in anchors]
                expected=[gf_mul(c,value) for c in matrix[i]]
                addr=coefficient_key(g,i)|value
                for plane in range(2):
                    k,l=other[plane*2:plane*2+2]
                    assert table_word(addr+plane*SLOTS,None)==(
                        displacement(expected[k])|(displacement(expected[l])<<32))
                for plane in range(2,4):
                    assert table_word(addr+plane*SLOTS,None)==sum(
                        expected[s]<<(TAG_BITS*(s%4)) for s in range(4*(plane-2),4*(plane-1)))
        for ranks in ((0,)*DIMENSION,(RADIUS-1,)*DIMENSION):
            for i in range(DIMENSION):
                assert table_word(descriptor_base(g,ranks)+i,None)==(
                    list_address(anchors[i],ranks[i])|((coefficient_key(g,i)+2*SLOTS)<<32))
    assert 48<=COEFF_BASE<COEFF_END<=DESCRIPTOR_BASE<DESCRIPTOR_END<=PARITY_BASE
    assert COEFF_BASE%ROW_STRIDE==0 and ROW_STRIDE==4*SLOTS
    assert PARITY_BASE+(1<<PRIVATE_BITS)==TAG_BASE
    assert TAG_BASE+SLOTS==LOW_SPLIT
    assert LOW_SPLIT+(1<<16)==HIGH_SPLIT and HIGH_SPLIT+(1<<15)==LIST_BASE
    assert LIST_BASE+(1<<PRIVATE_BITS)*RADIUS==TAG_SCRATCH
    assert TAG_SCRATCH+SIDES*RADIUS==TAG_SCRATCH_END<=LEADER_BASE
    assert LEADER_BASE+(1<<PARITY_BITS)==1<<32
    basis,_=split_basis()
    assert [split_syndrome(c) for c in basis]==[1<<i for i in range(PARITY_BITS)]
    assert [linear_map(split_syndrome(1<<i),basis) for i in range(PARITY_BITS)]==[1<<i for i in range(PARITY_BITS)]
    assert len({displacement(x) for x in range(SLOTS)})==SLOTS
    regs=set()
    ops=Counter()
    for op,dest,*operands in schedule():
        regs.update(x for x in (dest,*operands) if isinstance(x,str))
        if op!='move':
            ops[op]+=1
    W=sum(ops.values())
    assert W==EXPECTED_W and len(regs)<=64
    assert W**exponent.denominator<=1<<exponent.numerator
    return {"W":W,"L":str(exponent),"registers":len(regs),"ops":dict(ops),"preload_end":1<<32}


def main():
    from multiprocessing import get_context
    sector=json.load(sys.stdin)["sector"]
    assert 0<=sector<4
    data,history=rank_witness()
    rank_job=get_context('fork').Process(target=check_rank_sector,args=(sector,data[sector]))
    rank_job.start()
    hists=read_histograms()
    check_sector(sector,hists)
    profiles=profile_witness()
    kernels=local_pair_shells(hists,history)
    for changes in PROFILE_SECTORS[sector]:
        assert profile(changes,kernels,history)==profiles[changes]
    out={"sector":sector,"tag_rows":SECTOR_BOUNDARIES[sector+1]-SECTOR_BOUNDARIES[sector]}
    if sector==3:
        manifest=json.loads((HERE/'submission.json').read_text())
        B,L=[Fraction(int(c['score']['numerator']),int(c['score']['denominator'])) for c in manifest['claims']]
        n,occupied=component_numerator(profiles,history)
        assert 1-Fraction(n**COMPONENTS,10**2000)<=B
        out.update(check_program(L))
        out.update({"B":str(B),"last_occupied":str(occupied[-1]),"component_numerator":str(n)})
    rank_job.join()
    assert rank_job.exitcode==0
    out["exact_rank_weights"]=len(history)
    out["checked_profiles"]=PROFILE_SECTORS[sector]
    print(json.dumps(out))


if __name__=='__main__':
    main()
