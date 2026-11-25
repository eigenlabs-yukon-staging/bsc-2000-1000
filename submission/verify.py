"""Reproduce the fixed local spectra and the exact line-code certificate."""
import base64
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys
import zlib

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    PAIR_BASE, PAIR_PLANE, PAIR_STRIDE, PAIR_END, PAIRS, table_word,
    BLOCKS, COLUMNS, PRIVATE_BITS, PARITY_BITS, SIDE_BITS, SLOTS, LIST_COUNT,
    LEADER_BASE, PARITY_BASE, TAG_BASE, GAMMA_BASE, TAG_BITS, RECORD_BITS,
    CORRECTION_SHIFT, LIST_SHIFT, directions, displacement,
    gamma, schedule, split_basis, split_syndrome, linear_map, pair_slots, tag_tuple, syndrome_columns,
)

CUTOFF = 25
EXACT_DEPTH = 11
RANK_CAP = 64
SECTOR_BOUNDARIES = (0, 1100, 2200, 3300, 4096)
DIGIT_BITS = 78
SPECTRUM_MODULUS = 1 << ((2 * CUTOFF + 1) * DIGIT_BITS)


def gray_order(columns):
    values = [0]
    for column in columns:
        values.extend(value ^ column for value in reversed(values))
    return values


def pack(sequence, size):
    return sum(coefficient << (size * degree)
               for degree, coefficient in enumerate(sequence))


def coefficient(value, index, size):
    return (value >> (index * size)) & ((1 << size) - 1)


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


def outer_spectra(hists):
    evaluation_scale = 1 << SIDE_BITS
    max_evaluation = max(evaluation_scale, max(
        sum(value << (SIDE_BITS-i) for i,value in enumerate(row)) for row in hists
    ))
    assert SLOTS**2 * (1 << (2*CUTOFF)) * max_evaluation**4 < (1 << DIGIT_BITS) * evaluation_scale**4
    shift=7
    assert all(not any(row[:shift]) for row in hists[1:])
    triple_degree=2*CUTOFF-3*shift
    quadruple_degree=2*CUTOFF-4*shift
    triple_mod=1 << ((triple_degree+1)*DIGIT_BITS)
    quad_mod=1 << ((quadruple_degree+1)*DIGIT_BITS)
    triple=[0]+[pack(row[shift:shift+triple_degree+1],DIGIT_BITS) for row in hists[1:]]
    quad=[n&(quad_mod-1) for n in triple]
    four=sum(pow(n,4,quad_mod) for n in quad[1:])&(quad_mod-1)
    three=0
    lines=subspaces=0
    for s,d in directions():
        subspaces+=1
        three=(three+12*triple[s]*triple[d]*triple[s^d])&(triple_mod-1)
        for a in range(1,SLOTS):
            if a>min(a^s,a^d,a^s^d):
                continue
            lines+=1
            first=(quad[a]*quad[a^s])&(quad_mod-1)
            second=(quad[a^d]*quad[a^s^d])&(quad_mod-1)
            four=(four+12*first*second)&(quad_mod-1)
    assert subspaces==(SLOTS-1)//3
    assert lines==subspaces*(SLOTS//4-1)
    assert 1+(SLOTS-1)+12*(lines+subspaces)==SLOTS**2
    return [
        [0]*(r*shift)+[coefficient(n,w,DIGIT_BITS) for w in range(2*CUTOFF+1-r*shift)]
        for r,n in ((4,four),(3,three))
    ]


def cumulative_binomial(width):
    prefix = []
    for n in range(width + 1):
        row = []
        for w in range(n + 1):
            row.append((row[-1] if row else 0) + comb(n, w))
        prefix.append(row)
    return prefix


def pair_bounds(width, spectrum):
    prefix = cumulative_binomial(width)
    totals = []
    for t in range(CUTOFF + 1):
        pairs = 0
        for d in range(1, 2 * t + 1):
            if spectrum[d]:
                common = 0
                for a in range(max(0, d - t), min(d, t) + 1):
                    rest = min(width - d, t - max(a, d - a))
                    common += comb(d, a) * prefix[width - d][rest]
                pairs += common * spectrum[d]
        totals.append(pairs)
    return totals


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


def first_shell_pairs(weight, distance):
    lo = max((distance+1)//2, weight-(SIDE_BITS-distance))
    hi = min(distance, weight)
    return sum(comb(distance,a) * comb(SIDE_BITS-distance,weight-a)
               for a in range(lo,hi+1))


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


def shortlisted_volumes(hists, history):
    local=coverage(cumulative_binomial(SIDE_BITS)[SIDE_BITS][:CUTOFF+1],
                   pair_bounds(SIDE_BITS,hists[0]))
    assert local[-1]<=1<<31
    tagged=[sum(row[d] for row in hists[1:]) for d in range(2*CUTOFF+1)]
    combined=pair_bounds(SIDE_BITS,tagged)
    total=0
    rejected=0
    accept_lower=0
    accepted=0
    A=[]
    U=[]
    for w in range(CUTOFF+1):
        old_total=total
        if w<=EXACT_DEPTH:
            row=history[w]
            total+=sum(row)
            rejected+=sum(row[LIST_COUNT:])
        else:
            total=max(total,local[w])
            shell_pairs=sum(tagged[d]*first_shell_pairs(w,d)
                            for d in range(1,2*w+1))
            rejected=min(rejected+min(comb(SIDE_BITS,w),shell_pairs//LIST_COUNT),
                         combined[w]//(4*LIST_COUNT-2))
        accept_lower=max(accept_lower,total-rejected)
        nxt=min(accepted+total-old_total,accept_lower)
        A.append(nxt-accepted)
        U.append(total-old_total-A[-1])
        accepted=nxt
    assert all(x>=0 for x in A+U) and A[0]==1
    p4=convolution([A,A,A,A])
    p3=convolution([A,A,A,U])
    p2=convolution([A,A,U,U])
    volume=0
    output=[]
    for a,b,c in zip(p4,p3,p2,strict=True):
        volume+=a+4*b+6*c
        output.append(volume)
    return output


def component_numerator(hists, history):
    volumes=shortlisted_volumes(hists,history)
    four,three=outer_spectra(hists)
    quad_pairs=pair_bounds(4*SIDE_BITS,four)
    triple_pairs=pair_bounds(3*SIDE_BITS,three)
    pairs=[value+sum(min(comb(SIDE_BITS,w),1<<31)*triple_pairs[t-w]
                    for w in range(t+1))
           for t,value in enumerate(quad_pairs)]
    occupied=coverage(volumes,pairs)
    assert occupied[-1] <= 1 << 100
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
    assert TAG_BASE+SLOTS==GAMMA_BASE<=LEADER_BASE
    assert LIST_SHIFT==8 and RECORD_BITS==TAG_BITS+6
    assert LIST_SHIFT + RECORD_BITS*LIST_COUNT <= CORRECTION_SHIFT
    assert CORRECTION_SHIFT+PRIVATE_BITS==64
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
    assert count == 4944 and len(registers)<=64
    assert count**exponent.denominator <= 1<<exponent.numerator
    return {"W":count,"ops":dict(ops),"registers":len(registers),"working_end":GAMMA_BASE,
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
    if sector==3:
        data=json.loads((HERE/"submission.json").read_text())
        B,L=[Fraction(int(x["score"]["numerator"]),int(x["score"]["denominator"]))
             for x in data["claims"]]
        numerator,occupied=component_numerator(hists,history)
        assert 1-Fraction(numerator**BLOCKS,10**2000)<=B
        report.update(check_program(L))
        report.update({"component_numerator":str(numerator),"final_occupied":str(occupied[-1]),"B":str(B),"L":str(L)})
    rank_check.join()
    assert rank_check.exitcode==0
    report["exact_rank_weights"]=len(history)
    print(json.dumps(report))


if __name__=="__main__":
    main()
