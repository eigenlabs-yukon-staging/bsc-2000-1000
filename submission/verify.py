"""Reproduce the fixed local spectra and the exact line-code certificate."""
import base64
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import random
import sys
import zlib

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    SIDE_LIST_COUNTS,
    PAIR_BASE, PAIR_PLANE, PAIR_STRIDE, PAIR_END, PAIRS, table_word,
    BLOCKS, COLUMNS, PRIVATE_BITS, PARITY_BITS, SIDE_BITS, SLOTS, LIST_COUNT,
    END, LEADER_BASE, PARITY_BASE, TAG_BASE, GAMMA_BASE, LOW_SPLIT, HIGH_SPLIT,
    LIST_BASE, LIST_END, TAG_SCRATCH, TAGS_IN_REGISTERS, TAG_BITS, directions, displacement,
    gamma, schedule, split_basis, split_syndrome, linear_map, pair_slots, tag_tuple, syndrome_columns,
)

CUTOFF = 25
EXACT_DEPTH = 11
RANK_CAP = 64
SECTOR_BOUNDARIES = (0, 1100, 2200, 3300, 4096)
SIDES = 4
DIMENSION = 2
GLOBAL_CUTOFF = CUTOFF
SCALE = 1<<24
STRIDE = 256
DIGIT_MASK = (1<<STRIDE)-1
POLY_MODULUS = 1<<(STRIDE*(2*GLOBAL_CUTOFF+1))


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


def convolution(rows):
    result = [1] + [0] * CUTOFF
    for row in rows:
        result = [sum(result[d] * row[t-d] for d in range(t+1))
                  for t in range(CUTOFF+1)]
    return result


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
    assert len(data)==4 and LIST_COUNT<=RANK_CAP
    assert all(len(sector)==EXACT_DEPTH+1 for sector in data)
    assert all(len(row)==RANK_CAP+1 and all(type(x) is int and x>=0 for x in row)
               for sector in data for row in sector)
    total=[[sum(counts) for counts in zip(*rows,strict=True)]
           for rows in zip(*data,strict=True)]
    return data,total


def check_rank_sector(sector, witness):
    assert exact_rank_sector(sector)==witness


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


def ceildiv(a,b):
    return -(-a//b)


def checked_pack(counts):
    assert all(0<=x<=DIGIT_MASK for x in counts)
    return sum(x<<(STRIDE*i) for i,x in enumerate(counts))


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
            cap=threshold
            low_rows=range(min(w,largest)+1)
            last_columns=[min(w-u,largest) for u in low_rows]
            for k in range(-1,largest+1):
                first_cost=marginal(k)
                if first_cost>=cap:
                    break
                residual=0
                for u,last in zip(low_rows,last_columns,strict=True):
                    if last>k:
                        row=bounds[u]
                        available=row[last]-(row[k] if k>=0 else 0)
                        residual+=min(totals[u],available)
                cap=min(cap,first_cost+residual)
            assert cap>=previous
            shells.append(cap-previous)
            previous=cap
        result.append(shells)
    return result


def summed_profile(changes, kernels, history):
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


def row_kernels(hists, totals):
    width = len(totals)
    factors = {}
    for u in range(width):
        for z in range(u, width):
            choices = []
            for d in range(max(1, z-u), u+z+1):
                a = (d+u-z)//2
                b = (u+z-d)//2
                if (d+u-z)%2 == 0 and 0<=a<=d and 0<=b<=SIDE_BITS-d:
                    choices.append((d, comb(d,a)*comb(SIDE_BITS-d,b)))
            factors[u,z]=choices
    result=[]
    for hist in hists[1:]:
        ambient = [[0]*width for _ in totals]
        for (u,z), choices in factors.items():
            ambient[u][z] = ambient[z][u] = sum(hist[d]*q for d,q in choices)
        out = []
        for u,row in enumerate(ambient):
            remaining = totals[u]
            for z,cap in enumerate(row):
                take = min(remaining,cap,totals[z])
                remaining -= take
                out.append(take)
        result.append(out)
    return result


def square_profile(j, kernels, totals):
    width = len(totals)
    freedom = j-DIMENSION
    scale = 256
    bounds = []
    for position in range(width**2):
        total = sum(root_ceil(scale**freedom * row[position]**j, freedom)
                    for row in kernels if row[position])
        bounds.append(root_ceil(scale**(j-freedom) * total**freedom, j))
    coefficient_bound = sum(bounds)**j * sum(totals)**(SIDES-j)
    byte_width = (coefficient_bound.bit_length()+7)//8
    assert byte_width>0
    y_stride = SIDES*(width-1)+1
    rows = [bounds[i*width:(i+1)*width] for i in range(width)]
    def pack(rows):
        return int.from_bytes(b''.join(b''.join(x.to_bytes(byte_width,'little') for x in row)
                        + bytes(byte_width*(y_stride-len(row))) for row in rows),'little')
    pair = pack(rows)
    common = pack([[0]*i+[x] for i,x in enumerate(totals)])
    c = GLOBAL_CUTOFF
    modulus = 1<<(byte_width*8*y_stride*(c+1))
    merged = (pow(pair,j,modulus)*pow(common,SIDES-j,modulus))%modulus
    packed = merged.to_bytes(byte_width*y_stride*(c+1),'little')
    prefix = []
    running=0
    def coeff(i,k):
        offset = byte_width*(i*y_stride+k)
        return int.from_bytes(packed[offset:offset+byte_width],'little')
    for t in range(c+1):
        running += coeff(t,t)+sum(coeff(i,t)+coeff(t,i) for i in range(t))
        prefix.append(ceildiv(running,scale**j))
    return prefix



def exact_three_row_profile(kernels, totals, directions):
    width=len(totals)
    y_stride=4*(width-1)+1
    whole_bound=3*len(directions)*max(map(sum,kernels))**3*sum(totals)
    byte_width=(whole_bound.bit_length()+7)//8

    def pack(rows):
        return int.from_bytes(b"".join(
            b"".join(x.to_bytes(byte_width,"little") for x in row)+
            bytes(byte_width*(y_stride-len(row))) for row in rows
        ),"little")

    encoded=[pack([row[i:i+width] for i in range(0,width*width,width)]) for row in kernels]
    common=pack([[0]*i+[x] for i,x in enumerate(totals)])
    limit=byte_width*y_stride*(GLOBAL_CUTOFF+1)
    mask=(1<<(8*limit))-1
    product=0
    for first,second,third in directions:
        assert all(0<index<=len(kernels) for index in (first,second,third))
        a,b,c=(encoded[index-1] for index in (first,second,third))
        product=(product+((a*b)&mask)*c)&mask
    product=(3*product*common)&mask
    data=product.to_bytes(limit,"little")

    def coeff(i,k):
        offset=byte_width*(i*y_stride+k)
        return int.from_bytes(data[offset:offset+byte_width],"little")

    out=[]
    running=0
    for t in range(GLOBAL_CUTOFF+1):
        running+=coeff(t,t)+sum(coeff(i,t)+coeff(t,i) for i in range(t))
        out.append(running)
    return out



def profile_witness():
    data=json.loads((HERE/'profiles.json').read_text())
    assert len(data)==2
    assert all(len(row)==GLOBAL_CUTOFF+1 for row in data)
    assert all(isinstance(x,str) and x.isdecimal() for row in data for x in row)
    return {support:list(map(int,row)) for support,row in enumerate(data,3)}


def verified_profile(support,hists,history):
    totals=[sum(row) for row in history]
    assert len(totals)==EXACT_DEPTH+1 and sum(totals)==1<<PARITY_BITS
    diagonal=summed_profile(support,local_pair_shells(hists,history),history)
    rows=row_kernels(hists,totals)
    if support==3:
        triples=[(a,b,a^b) for a,b in directions()]
        assert len(triples)==1365
        assert set(t for triple in triples for t in triple)==set(range(1,SLOTS))
        assert all(gamma(a)==b and gamma(b)==c and gamma(c)==a for a,b,c in triples)
        square=exact_three_row_profile(rows,totals,triples)
    else:
        assert support==4
        square=square_profile(support,rows,totals)
    return [min(a,b) for a,b in zip(diagonal,square,strict=True)]


def component_numerator(profiles,history):
    totals=[sum(row) for row in history]
    assert sum(totals)==1<<PARITY_BITS
    admitted=[[sum(row[:r]) for row in history] for r in SIDE_LIST_COUNTS]
    assert all(sum(row)==r*(1<<PRIVATE_BITS) for row,r in zip(admitted,SIDE_LIST_COUNTS,strict=True))
    rejected=[[t-a for t,a in zip(totals,row,strict=True)] for row in admitted]
    padding=[0]*(CUTOFF+1-len(totals))
    admitted=[row+padding for row in admitted]
    rejected=[row+padding for row in rejected]
    shell=[0]*(CUTOFF+1)
    for mask in range(1<<SIDES):
        if mask.bit_count()<DIMENSION:
            continue
        poly=convolution([admitted[i] if (mask>>i)&1 else rejected[i] for i in range(SIDES)])
        shell=[x+y for x,y in zip(shell,poly,strict=True)]
    from itertools import accumulate
    volumes=list(accumulate(shell))
    pairs=[4*a+b for a,b in zip(profiles[3],profiles[4],strict=True)]
    occupied=coverage(volumes,pairs)
    assert occupied[-1]<=1<<100
    numerator=sum((x-(occupied[t-1] if t else 0))*9**(200-t)
                  for t,x in enumerate(occupied))
    return numerator,occupied


def check_program(exponent):
    from collections import Counter
    from itertools import combinations
    assert 48<=PAIR_BASE<PAIR_END<=PARITY_BASE
    assert PAIR_BASE%PAIR_STRIDE==0 and PAIR_PLANE==SLOTS*SLOTS
    assert PAIR_STRIDE==4*PAIR_PLANE and PAIR_END==PAIR_BASE+len(PAIRS)*PAIR_STRIDE
    assert PARITY_BASE + (1 << PRIVATE_BITS) == TAG_BASE
    assert TAG_BASE+SLOTS==GAMMA_BASE and GAMMA_BASE+SLOTS==LOW_SPLIT
    assert HIGH_SPLIT-LOW_SPLIT==1<<16 and LIST_BASE-HIGH_SPLIT==1<<15
    assert LIST_END-LIST_BASE==(1<<PRIVATE_BITS)*LIST_COUNT
    assert LIST_END==TAG_SCRATCH
    assert END==TAG_SCRATCH+(0 if TAGS_IN_REGISTERS else 4*LIST_COUNT) <= LEADER_BASE
    assert LEADER_BASE+(1<<PARITY_BITS)==1<<32
    assert len({displacement(x) for x in range(SLOTS)})==SLOTS
    basis,_=split_basis()
    assert len(basis)==PARITY_BITS and len(COLUMNS)==PRIVATE_BITS+12
    for i,row in enumerate(basis):
        assert split_syndrome(row)==1<<i
    for i in range(PARITY_BITS):
        coords=split_syndrome(1<<i)
        assert linear_map(coords,basis)==1<<i
    for t in range(SLOTS):
        image=gamma(t)
        assert gamma(image)==(t^image)
        assert t==0 or image not in (0,t)
    for i,j in combinations(range(4),2):
        k,l=pair_slots(i,j)
        assert len({i,j,k,l})==4
        for s in range(SLOTS):
            tags=tag_tuple(0,s)
            delta=gamma(tags[i]^tags[j])
            assert tags[k]==tags[i]^delta and tags[l]==tags[j]^delta
    for pair_id,(i,j) in enumerate(PAIRS):
        k,l=pair_slots(i,j)
        for a in (0,1,3,1023,4095):
            for d in (0,1,31,2048,4095):
                tags=tag_tuple(a,d)
                key=PAIR_BASE+pair_id*PAIR_STRIDE+tags[i]*SLOTS+tags[j]
                assert table_word(key,None)==displacement(tags[k]) | (displacement(tags[l])<<32)
                assert table_word(key+PAIR_PLANE,None)==0
                assert table_word(key+2*PAIR_PLANE,None)==sum(t<<(TAG_BITS*s) for s,t in enumerate(tags))
    ops=Counter()
    registers=set()
    for op,dest,*operands in schedule():
        if isinstance(dest,str):
            registers.add(dest)
        registers.update(value for value in operands if isinstance(value,str))
        if op!='move':
            ops[op]+=1
    count=sum(ops.values())
    assert count == 16104 and len(registers)<=64
    assert count**exponent.denominator <= 1<<exponent.numerator
    return {"W":count,"ops":dict(ops),"registers":len(registers),"working_end":END,
            "preload_end":1<<32,"list_count":LIST_COUNT}


def main():
    sector=json.load(sys.stdin)["sector"]
    assert 0<=sector<4
    hists=read_histograms()
    from multiprocessing import get_context
    data,history=rank_witness()
    rank_check=get_context('fork').Process(target=check_rank_sector,args=(sector,data[sector]))
    rank_check.start()
    check_sector(sector,hists)
    report={"sector":sector,"tag_rows":SECTOR_BOUNDARIES[sector+1]-SECTOR_BOUNDARIES[sector]}
    profiles=profile_witness()
    if sector in (1,2):
        support=sector+2
        assert verified_profile(support,hists,history)==profiles[support]
        report['checked_support']=support
    if sector==3:
        data=json.loads((HERE/"submission.json").read_text())
        B,L=[Fraction(int(x["score"]["numerator"]),int(x["score"]["denominator"]))
             for x in data["claims"]]
        numerator,occupied=component_numerator(profiles,history)
        assert 1-Fraction(numerator**BLOCKS,10**2000)<=B
        report.update(check_program(L))
        report.update({"component_numerator":str(numerator),"final_occupied":str(occupied[-1]),"B":str(B),"L":str(L)})
    rank_check.join()
    assert rank_check.exitcode==0
    report["exact_rank_weights"]=len(history)
    print(json.dumps(report))


if __name__=="__main__":
    main()
