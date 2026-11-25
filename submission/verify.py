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
from collections import Counter
from functools import lru_cache
from construction import (
    PRIVATE_BITS, PARITY_BITS, SIDE_BITS, TAG_BITS, SLOTS, COLUMNS,
    SIDES, DIMENSION, COMPONENT_BITS, PAYLOAD_BITS, COMPONENTS, RADIUS,
    MODULUS, EVALUATION_POINTS, ANCHORS, FIELD_DISP_BASE, FIELD_TAG_BASE,
    GROUP_BASE, GROUP_END, GROUP_STRIDE, PARITY_BASE, TAG_BASE, LOW_SPLIT,
    HIGH_SPLIT, LIST_BASE, LIST_END, LEADER_BASE, U_WORK, S_WORK, SHORTLIST,
    WEIGHT_WORK, DISP_WORK, PREFIX_WORK, SPILLED_DEPTHS, BEST_SCORE, BEST_KEY, WORK_END, PAIRED,
    gf_mul, gf_pow, coefficients, table_word, displacement, split_basis,
    split_syndrome, syndrome_columns, linear_map, allocate_registers,
    group_schedule, setup_component, decode_winner, choose_extra, make_priorities, PRIORITY_WORK, EXTRA_KEY, KEY_RADIX,
)

EXACT_DEPTH = 11
RANK_CAP = 64
SECTOR_BOUNDARIES = (0, 1100, 2200, 3350, 4096)


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
    populations=[]
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
        from array import array
        populations.append(array('H',previous))
        ranks=[]
        running=0
        for x in delta[:-1]:
            running+=x
            ranks.append(running)
        ranks.append(above)
        history.append(ranks)
    classes={}
    for q in range(size):
        w2=next(w for w,population in enumerate(populations) if population[q]>=2)
        w3=next(w for w,population in enumerate(populations) if population[q]>=3)
        w4=next(w for w,population in enumerate(populations) if population[q]>=4)
        priority=(w3<<10)|((31-w2)<<5)|(31-w4)
        counters=classes.setdefault(str(priority),[[0,0,0] for _ in populations])
        old=0
        for w,population in enumerate(populations):
            new=population[q]
            counters[w][0]+=min(new,2)-min(old,2)
            counters[w][1]+=int(old<3<=new)
            counters[w][2]+=new-old
            old=new
    return history,classes


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


def check_rank_sector(sector,witness,priority):
    assert exact_rank_sector(sector)==(witness,priority)


EXPECTED_W=30985433374
PROFILE_SECTORS=((11,12,13),(14,15,16),(17,18),(19,20))

GLOBAL_CUTOFF=135
SCALE=1<<24
STRIDE=1280
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


def profile(changes, kernels, history):
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


def residue_profile(j,kernels,totals,sides,dimension,cutoff,root_ceil,modulus,scale):
    width=len(totals)
    freedom=j-dimension
    exponent_span=2*sides*(width-1)+1
    changed=[[0]*exponent_span for _ in range(modulus)]
    unchanged=[[0]*exponent_span for _ in range(modulus)]
    for position in range(width**2):
        first,second=divmod(position,width)
        sum_weights=first+second
        total=sum(root_ceil(scale**freedom*row[position]**j,freedom)
                  for row in kernels if row[position])
        norm=root_ceil(scale**(j-freedom)*total**freedom,j)
        changed[second%modulus][sum_weights]+=-(-norm//(1<<sum_weights))
    for u,t in enumerate(totals):
        unchanged[u%modulus][2*u]=-(-(scale*t)//(1<<(2*u)))
    bound=sum(map(sum,changed))**j*sum(map(sum,unchanged))**(sides-j)
    byte_width=(bound.bit_length()+7)//8
    box_bytes=byte_width*exponent_span*modulus
    box_bits=8*box_bytes
    mask=(1<<box_bits)-1

    def pack(rows):
        return int.from_bytes(b''.join(x.to_bytes(byte_width,'little')
                                       for row in rows for x in row),'little')

    def multiply(a,b):
        merged=a*b
        while merged>>box_bits:
            merged=(merged&mask)+(merged>>box_bits)
        return merged

    def power(poly,count):
        result=1
        while count:
            if count&1:
                result=multiply(result,poly)
            count>>=1
            if count:
                poly=multiply(poly,poly)
        return result

    product=multiply(power(pack(changed),j),power(pack(unchanged),sides-j))
    data=product.to_bytes(box_bytes,'little')

    def coeff(residue,index):
        start=byte_width*(residue*exponent_span+index)
        return int.from_bytes(data[start:start+byte_width],'little')

    decoded=[[coeff(r,d)<<d for r in range(modulus)] for d in range(2*cutoff+1)]
    degree_sums=[sum(row) for row in decoded]
    result=[]
    for t in range(cutoff+1):
        count=0
        for d,row in enumerate(decoded[:2*t+1]):
            start=max(0,d-t)
            stop=min(t,d)+1
            if stop-start>=modulus:
                count+=degree_sums[d]
            else:
                count+=sum(row[v%modulus] for v in range(start,stop))
        result.append(-(-count//scale**sides))
    return result


def profile_witness():
    data=json.loads((HERE/'profiles.json').read_text())
    assert len(data)==SIDES-DIMENSION
    assert all(len(row)==GLOBAL_CUTOFF+1 for row in data)
    assert all(isinstance(x,str) and x.isdecimal() for row in data for x in row)
    return {changes:list(map(int,row))
            for changes,row in enumerate(data,DIMENSION+1)}


def component_numerator(profiles, history, priority):
    from decimal import Decimal, Inexact, Rounded, localcontext

    totals={}
    for sector in priority:
        for key,rows in sector.items():
            if key not in totals:
                totals[key]=[[0,0,0] for _ in history]
            for accumulated,row in zip(totals[key],rows,strict=True):
                for i,n in enumerate(row):
                    accumulated[i]+=n
    observed=[[sum(rows[w][i] for rows in totals.values()) for i in range(3)]
              for w in range(len(history))]
    assert observed==[[sum(row[:2]),row[2],sum(row)] for row in history]
    assert sum(row[2] for row in observed)==1<<PARITY_BITS

    # Pack a bivariate polynomial in top-two count and weight.
    degree=(len(history)-1)*SIDES+1
    digits=190
    assert (1<<PARITY_BITS)**SIDES<10**digits

    def encoded(rows,kind):
        counts={}
        for w,(a,b,total) in enumerate(rows):
            if kind=='plain':
                counts[w]=total-a
                counts[degree+w]=a
            elif kind=='disabled':
                counts[w]=total-a-b
            else:
                counts[w]=b
        return Decimal(''.join(str(counts.get(i,0)).zfill(digits)
                               for i in range(degree+len(rows)-1,-1,-1)))

    with localcontext() as context:
        context.prec=digits*(SIDES+1)*degree+20
        context.Emax=9999999
        context.traps[Inexact]=True
        context.traps[Rounded]=True
        plain={key:encoded(rows,'plain') for key,rows in totals.items()}
        disabled={key:encoded(rows,'disabled') for key,rows in totals.items()}
        third={key:encoded(rows,'third') for key,rows in totals.items()}
        all_plain=sum(plain.values())
        base=format(all_plain**SIDES,'f')
        later=all_plain
        lower=Decimal(0)
        added=Decimal(0)
        for key in sorted(totals,key=int):
            later-=plain[key]
            before=lower+disabled[key]+later
            after=lower+plain[key]+later
            divisor=after-before
            difference=after**SIDES-before**SIDES
            assert divisor>0 and third[key]>0
            quotient=difference//divisor
            assert quotient*divisor==difference
            added+=third[key]*quotient
            lower+=disabled[key]
        extension=format(added,'f')

    def coefficient_at(poly,index):
        stop=len(poly)-digits*index
        return int(poly[max(0,stop-digits):stop]) if stop>0 else 0

    assert sum(coefficient_at(base,i) for i in range((SIDES+1)*degree))==(1<<PARITY_BITS)**SIDES
    cumulative=0
    volume=[]
    for w in range(GLOBAL_CUTOFF+1):
        cumulative+=sum(coefficient_at(base,k*degree+w) for k in range(DIMENSION,SIDES+1))
        cumulative+=coefficient_at(extension,(DIMENSION-1)*degree+w)
        volume.append(cumulative)
    pairs=[sum(comb(SIDES,c)*profiles[c][t] for c in profiles)
           for t in range(GLOBAL_CUTOFF+1)]
    occupied=coverage(volume,pairs)
    assert occupied[-1]<=1<<(COMPONENT_BITS-PAYLOAD_BITS)
    p=sum((x-(occupied[t-1] if t else 0))*9**(COMPONENT_BITS-t)
          for t,x in enumerate(occupied))
    return p,occupied


def polynomial_gcd(a,b):
    while b:
        while a.bit_length()>=b.bit_length():
            a^=b<<(a.bit_length()-b.bit_length())
        a,b=b,a
    return a


def count_segment(instructions):
    ops=0
    regs=set()
    for op,dest,*args in allocate_registers(instructions):
        ops+=op!='move'
        regs.update(x for x in (dest,*args) if isinstance(x,str))
    return ops,len(regs)


def count_work():
    radices=[2]*(DIMENSION-1)+[3]
    precompute=4*PAIRED+sum(radices)*(4+7*PAIRED)
    nodes=1
    tree=0
    for depth,radix in enumerate(radices):
        nodes*=radix
        ops=PAIRED+1+PAIRED+int(depth>0)
        if 0<depth<=SPILLED_DEPTHS:
            ops+=PAIRED+1
        if depth<SPILLED_DEPTHS:
            ops+=PAIRED+1
        tree+=nodes*ops
    ordinary=precompute+tree+nodes*(PAIRED*6+1+3)+4
    registers=[]
    for g in (0,len(ANCHORS)-1):
        anchors=ANCHORS[g]
        for i in range(DIMENSION):
            observed=count_segment(group_schedule(g,anchors,i))
            assert observed[0]==ordinary-(0 if g else 5)
            registers.append(observed[1])
        branch_select=count_segment(choose_extra(anchors))
        assert branch_select[0]==38
        registers.append(branch_select[1])
    prep=count_segment(make_priorities())
    assert prep[0]==260
    registers.append(prep[1])
    limit=ordinary+38+2+2*(DIMENSION-1)
    offsets=[(count_segment(setup_component(b)),count_segment(decode_winner(b)))
             for b in range(COMPONENTS)]
    W=sum(a[0]+b[0]+prep[0]+len(ANCHORS)*limit-5 for a,b in offsets)
    used=max(registers+[p[1] for pair in offsets for p in pair])
    return W,used,offsets,ordinary,limit


def check_program(exponent):
    assert SIDES*EXACT_DEPTH < 256
    assert MODULUS.bit_length()==TAG_BITS+1
    assert gf_pow(2,SLOTS)==2
    for p in (2,3):
        assert polynomial_gcd(MODULUS,gf_pow(2,1<<(TAG_BITS//p))^2)==1
    assert len(set(EVALUATION_POINTS))==SIDES
    assert COMPONENT_BITS*COMPONENTS==2000 and PAYLOAD_BITS*COMPONENTS==1000
    assert len(ANCHORS)==comb(SIDES,DIMENSION)==184756
    assert RADIUS==4 and KEY_RADIX==3 and SIDES==2*DIMENSION and PAIRED*2==DIMENSION
    assert 48<=FIELD_DISP_BASE<FIELD_TAG_BASE<GROUP_BASE<GROUP_END<=PARITY_BASE
    assert FIELD_TAG_BASE-FIELD_DISP_BASE==SLOTS*SLOTS
    assert FIELD_TAG_BASE+SLOTS*SLOTS<=GROUP_BASE
    assert GROUP_STRIDE==DIMENSION*(SIDES+1)
    assert GROUP_END-GROUP_BASE==len(ANCHORS)*GROUP_STRIDE
    assert list((U_WORK,S_WORK,SHORTLIST,WEIGHT_WORK,DISP_WORK,PREFIX_WORK,BEST_SCORE,BEST_KEY,PRIORITY_WORK,EXTRA_KEY,WORK_END))==[
        LIST_END,LIST_END+SIDES,LIST_END+2*SIDES,LIST_END+2*SIDES+SIDES*RADIUS,
        LIST_END+2*SIDES+SIDES*RADIUS+DIMENSION*RADIUS,
        LIST_END+2*SIDES+SIDES*RADIUS+DIMENSION*RADIUS*(PAIRED+1),
        PREFIX_WORK+SPILLED_DEPTHS*(PAIRED+1),
        BEST_SCORE+1,BEST_KEY+1,PRIORITY_WORK+SIDES,EXTRA_KEY+1,
    ]
    assert PARITY_BASE+(1<<PRIVATE_BITS)==TAG_BASE and TAG_BASE+SLOTS==LOW_SPLIT
    assert LOW_SPLIT+(1<<16)==HIGH_SPLIT and HIGH_SPLIT+(1<<15)==LIST_BASE
    assert LIST_BASE+(1<<PRIVATE_BITS)*RADIUS==LIST_END
    assert WORK_END<=LEADER_BASE and LEADER_BASE+(1<<PARITY_BITS)==1<<32
    basis,_=split_basis()
    assert [split_syndrome(c) for c in basis]==[1<<i for i in range(PARITY_BITS)]
    assert [linear_map(split_syndrome(1<<i),basis) for i in range(PARITY_BITS)]==[1<<i for i in range(PARITY_BITS)]
    assert len({displacement(x) for x in range(SLOTS)})==SLOTS

    for g in (0,1,len(ANCHORS)//2,len(ANCHORS)-1):
        matrix=coefficients(ANCHORS[g])
        for i in range(DIMENSION):
            assert tuple(matrix[i][q] for q in ANCHORS[g])==tuple(int(i==j) for j in range(DIMENSION))
            base=GROUP_BASE+GROUP_STRIDE*g+(SIDES+1)*i
            assert table_word(base,None)==SHORTLIST+RADIUS*ANCHORS[g][i]
            for q,c in enumerate(matrix[i]):
                assert table_word(base+1+q,None)==FIELD_TAG_BASE+c*SLOTS
                for value in (0,1,301,4095):
                    x=gf_mul(c,value)
                    assert table_word(FIELD_TAG_BASE+c*SLOTS+value,None)==x
                    assert table_word(FIELD_DISP_BASE+c*SLOTS+value,None)==displacement(x)
    W,regs,offsets,ordinary,limit=count_work()
    assert W==EXPECTED_W and regs<=64
    assert W**exponent.denominator<=1<<exponent.numerator
    return {"W":W,"L":str(exponent),"registers":regs,"stage_endpoints":offsets,
            "ordinary_case":ordinary,"group_work_bound":limit,"preload_end":1<<32}


def main():
    from multiprocessing import get_context
    sector=json.load(sys.stdin)["sector"]
    assert 0<=sector<4
    data,history=rank_witness()
    priority=json.loads((HERE/"priority.json").read_text())
    assert len(priority)==4
    rank_job=get_context('fork').Process(target=check_rank_sector,args=(sector,data[sector],priority[sector]))
    rank_job.start()
    hists=read_histograms()
    check_sector(sector,hists)
    profiles=profile_witness()
    kernels=local_pair_shells(hists,history)
    totals=[sum(row) for row in history]
    rows=row_kernels(hists,totals)
    for changes in PROFILE_SECTORS[sector]:
        diagonal=profile(changes,kernels,history)
        residue=residue_profile(changes,rows,totals,SIDES,DIMENSION,GLOBAL_CUTOFF,root_ceil,16,4096)
        assert [min(a,b) for a,b in zip(diagonal,residue,strict=True)]==profiles[changes]
    out={"sector":sector,"tag_rows":SECTOR_BOUNDARIES[sector+1]-SECTOR_BOUNDARIES[sector]}
    if sector==3:
        manifest=json.loads((HERE/'submission.json').read_text())
        B,L=[Fraction(int(c['score']['numerator']),int(c['score']['denominator'])) for c in manifest['claims']]
        n,occupied=component_numerator(profiles,history,priority)
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
