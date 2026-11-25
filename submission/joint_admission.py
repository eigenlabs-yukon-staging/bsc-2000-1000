"""Lower dyadic coefficients for true-path admission by exact noise weight."""
from math import comb

from construction import ANCHOR_BITS, BUDGETS, LIMIT, M, cost
from dynamics import SCHEDULE, data, derive, kept_counts

PROB_BITS = 48
PROB_SCALE = 1<<PROB_BITS
STRIDE = 104
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


def child_removed_counts(cutoff):
    child_bits=2000-ANCHOR_BITS-M
    channel_denominator=10**child_bits
    unit=[PROB_SCALE*9**(child_bits-w) for w in range(cutoff+1)]
    all_children=[comb(child_bits,w) for w in range(cutoff+1)]
    local=[0]*(cutoff+1)
    for d,polynomial in enumerate(suffix_lower(cutoff)):
        anchors=comb(ANCHOR_BITS,d)
        for w in range(cutoff-d+1):
            low=polynomial&COEFFICIENT_MASK
            admitted=(low*channel_denominator+unit[w]-1)//unit[w]
            assert 0<=admitted<=all_children[w]
            local[d+w]+=anchors*(all_children[w]-admitted)
            polynomial>>=STRIDE
    remaining=[comb(M,w) for w in range(cutoff+1)]
    return [sum(value*remaining[t-s] for s,value in enumerate(local[:t+1]))
            for t in range(cutoff+1)]
