"""Lower dyadic coefficients for true-path admission by exact noise weight."""
from fractions import Fraction
from math import comb

from construction import ANCHOR_BITS, BUDGETS, GROUPS, LIMIT, M, cost
from dynamics import SCHEDULE, data, derive, kept_counts

PROB_BITS = 48
PROB_SCALE = 1<<PROB_BITS
STRIDE = 104
OUTPUT_SCALE = 10**60
COEFFICIENT_MASK = (1<<(PROB_BITS+1))-1


def local_weights():
    result = {}
    for width,histogram in data():
        noise = [9**(2*width-w) for w in range(2*width+1)]
        denominator = 10**(2*width)
        rows = []
        for profile,multiplicity in histogram.items():
            normalizer = sum(n*q for n,q in zip(profile,noise,strict=True))
            terms = []
            for w,(kept,q) in enumerate(zip(kept_counts(width,profile),noise,strict=True)):
                mass = (multiplicity*kept*q*PROB_SCALE)//denominator
                if mass:
                    terms.append((cost(normalizer,q),STRIDE*w,mass))
            rows.append(terms)
        result[width] = rows
    result[0] = [[(cost(10,9),0,9*PROB_SCALE//10),
                  (cost(10,1),STRIDE,PROB_SCALE//10)]]
    return result


def suffix_lower(cutoff):
    mask = (1<<(STRIDE*(cutoff+1)))-1
    lanes = sum(COEFFICIENT_MASK<<(STRIDE*w) for w in range(cutoff+1))
    types = local_weights()
    _,caps,_ = derive(10,True)
    state_count = LIMIT+1
    suffix = [PROB_SCALE]*state_count
    for index in range(len(SCHEDULE)-1,-1,-1):
        width = SCHEDULE[index]
        reserve = len(SCHEDULE)-1-index
        numerator = [0]*state_count
        for terms,cap_row in zip(types[width],caps[index],strict=True):
            by_cap = [0]*state_count
            for price,shift,mass in terms:
                for cap in range(price+reserve,state_count):
                    by_cap[cap] += (suffix[cap-price]*mass)<<shift
            for r,cap in enumerate(cap_row):
                numerator[r] += by_cap[cap]
        suffix = [((value&mask)>>PROB_BITS)&lanes for value in numerator]
    return [suffix[cap] for cap in BUDGETS]


def child_miss_upper(alpha):
    cutoff = len(alpha)-1
    child_bits = 2*GROUPS*M-ANCHOR_BITS-M
    child_den = 10**child_bits
    exterior_den = 10**M
    exterior = [comb(M,w)*9**(M-w) for w in range(min(M,cutoff)+1)]
    ext_weight = []
    for s in range(cutoff+1):
        weighted = sum(q*a for q,a in zip(exterior,alpha[s:]))
        scaled = weighted*OUTPUT_SCALE/exterior_den
        ext_weight.append(-(-scaled.numerator//scaled.denominator))
    child_weight = [comb(child_bits,z)*9**(child_bits-z)*PROB_SCALE
                    for z in range(cutoff+1)]
    numerator = 0
    for d,poly in enumerate(suffix_lower(cutoff)):
        anchor = comb(ANCHOR_BITS,d)*9**(ANCHOR_BITS-d)
        for z in range(cutoff-d+1):
            admitted = poly&COEFFICIENT_MASK
            miss = child_weight[z]-admitted*child_den
            assert miss>=0
            numerator += anchor*miss*ext_weight[d+z]
            poly >>= STRIDE
    denominator = 10**ANCHOR_BITS*child_den*PROB_SCALE
    return (numerator+denominator-1)//denominator
