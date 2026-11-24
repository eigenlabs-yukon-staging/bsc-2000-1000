"""Exact double-direction occupancy, admitted shells, and single collisions."""
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    AFFINE_BASE, AGREE, BANDS, CAP_BASE, DEPTH_BASE, DISAGREE, MASK, M, MEMORY_END, MODULUS,
    MULT, RADIUS, REDUCE_4, SKIP, SLICE, SLICES, SQUARE_BASE, affine_entry, cube, encode, inner_cap,
    operation_bound, poly_remainder, reduction_entry, square, square_entry, table_word,
)

CUTOFF = 230
Q = 1 << (2 * M)
FIELD_ORDER = 1 << M
COUNTS = tuple(comb(M, i) for i in range(CUTOFF + 1))


@lru_cache(None)
def binomial_prefix(n):
    row = [0]
    for k in range(n+1):
        row.append(row[-1] + comb(n, k))
    return row


def sorted_groups(values):
    groups = sorted(values, reverse=True)
    limits = [0]
    masses = [0]
    for value, count in groups:
        limits.append(limits[-1] + count)
        masses.append(masses[-1] + value * count)
    return groups, limits, masses


def leading_bias(table, capacity):
    groups, limits, masses = table
    i = bisect_left(limits, capacity) - 1
    return masses[i] + (capacity - limits[i]) * groups[i][0]


def double_direction_bounds():
    signed = [[] for _ in range(CUTOFF + 1)]
    for weight in range(1, 2*M + 1):
        gap = 2*M - 1 - 2*(weight-1)
        prev, value = 0, 1
        one = 3*comb(M, weight)
        two = 9*(comb(2*M, weight) - 2*comb(M, weight))
        assert one >= 0 and two >= 0
        for radius, row in enumerate(signed):
            row.append((abs(value), one, two))
            if radius < CUTOFF:
                n, rem = divmod(gap*value - (2*M-radius)*prev, radius+1)
                assert rem == 0
                prev, value = value, n
    one = [sorted_groups((v,c) for v,c,_ in row if c) for row in signed]
    two = [sorted_groups((v,c) for v,_,c in row if c) for row in signed]
    one_bias = [[leading_bias(table, c) for c in COUNTS] for table in one]
    ball = binomial_prefix(2*M)
    result = []
    for t in range(CUTOFF+1):
        disparity = 0
        for a in range(min(t,RADIUS)+1):
            na = COUNTS[a]
            for b in range(a, t-a+1):
                nb = COUNTS[b]
                radius = t-a-b
                x = na*one_bias[radius][b] + nb*one_bias[radius][a]
                x += leading_bias(two[radius], na*nb) + na*nb*ball[radius+1]
                disparity += x if a==b else 2*x
        result.append(disparity)
    assert all(a <= b for a,b in zip(result,result[1:]))
    return result


def single_direction_weights():
    coefficients = [
        sum(comb(M, v)*COUNTS[w-2*v] for v in range(w//2+1))
        for w in range(CUTOFF+1)
    ]
    weights = [0]*(CUTOFF+1)
    for j in range(CUTOFF+1):
        volume = 0
        for w in range(CUTOFF+1-j):
            count = coefficients[w] << j
            volume += count
            weights[j+w] += COUNTS[j]*count*min(FIELD_ORDER,2*volume)
        for w in range(2,len(coefficients)):
            coefficients[w] -= coefficients[w-2]
            assert coefficients[w] >= 0
    return convolution(weights,COUNTS)


@lru_cache(None)
def outer_counts(low, high):
    return [
        sum(COUNTS[i]*COUNTS[w-i] for i in range(low,min(high,w)+1))
        for w in range(CUTOFF+1)
    ]


@lru_cache(None)
def inner_counts(budget, low, high):
    result = [0]*(CUTOFF+1)
    for j in range(min(CUTOFF,budget)+1):
        cap = inner_cap(j,budget)
        prefix = binomial_prefix(j)
        for i in range(min(cap,(CUTOFF-j)//2)+1):
            until = min(j+1,max(0, high-i+1))
            before = min(j+1,max(0, low-i))
            choices = prefix[until]-prefix[before]
            result[j+2*i] += COUNTS[j]*comb(M-j,i)*choices
    return result


def convolution(a, b):
    return [sum(a[i]*b[t-i] for i in range(t+1)) for t in range(CUTOFF+1)]


def admitted_shells():
    either_orientation = [0]*(CUTOFF+1)
    both = [0]*(CUTOFF+1)
    for i,(lo,hi,budget) in enumerate(BANDS):
        prod = convolution(outer_counts(lo,hi), inner_counts(budget,0,M))
        either_orientation = [x+y for x,y in zip(either_orientation,prod,strict=True)]
        for j in range(i,len(BANDS)):
            other_lo,other_hi,other_budget = BANDS[j]
            prod = convolution(inner_counts(other_budget,lo,hi),
                               inner_counts(budget,other_lo,other_hi))
            multiplier = 1 if i == j else 2
            both = [x+multiplier*y for x,y in zip(both,prod,strict=True)]
    result = [2*x-y for x,y in zip(either_orientation,both,strict=True)]
    for t,count in enumerate(result):
        assert 0 <= count <= comb(4*M,t)
    return result


def certificate():
    pairs = double_direction_bounds()
    shells = admitted_shells()
    lower = []
    previous = volume = success = 0
    penalty = Fraction(0)
    weighted_single = single_direction_weights()
    for t,(count,collision) in enumerate(zip(shells,pairs,strict=True)):
        volume += count
        k = 1 + collision//Q
        alpha = Fraction(2,k+1)
        current = max(previous, volume*(2*k*Q-collision)//(Q*k*(k+1)))
        assert current <= Q
        wt = 9**(4*M-t)
        success += (current-previous)*wt
        penalty += weighted_single[t]*wt*alpha
        lower.append(current)
        previous = current
    removal = 2*penalty/(FIELD_ORDER*10**(4*M))
    return 1-Fraction(success,10**(4*M))+removal,success,removal,lower


def field_checks():
    x = 2
    checks = {}
    for i in range(1,M+1):
        x = square(x)
        if i in (100,250):
            a,b = x^2,MODULUS
            while b:
                a,b=b,poly_remainder(a,b)
            assert a == 1
            checks[i] = a
    assert x == 2
    rng = random.Random(93451)
    for x in [0,1,MASK] + [rng.getrandbits(M) for _ in range(10)]:
        y = 0
        for j in range((M+7)//8):
            y ^= square_entry(j,(x>>(8*j))&255)
        assert y == square(x)
    assert max(reduction_entry(x) for x in range(16)).bit_length() <= 64
    assert SLICES*SLICE <= 512
    for i in (0,13,14,63,64,127,499):
        x = rng.getrandbits(M)
        delta = 0
        for j in range(SLICES):
            delta ^= affine_entry(i,j,(x>>(SLICE*j))&((1<<SLICE)-1))
        assert cube(x)^delta == cube(x^(1<<i))
    return checks


def table_checks():
    assert BANDS[0][0] == 0
    for left,right in zip(BANDS,BANDS[1:]):
        assert left[1]+1 == right[0]
    assert 0 <= RADIUS <= M
    assert 48 < DISAGREE < AGREE < MULT < REDUCE_4 < SQUARE_BASE < DEPTH_BASE < CAP_BASE < AFFINE_BASE < MEMORY_END <= 1 << 32
    for lo,hi,power in BANDS:
        base = table_word(DEPTH_BASE+lo)
        assert table_word(DEPTH_BASE+hi) == base
        for j in range(M+1):
            cap = table_word(base+j)
            if cap == SKIP:
                assert j > power
            else:
                n = M-j
                actual = binomial_prefix(n)[cap+1]
                assert actual << j <= 1 << power
                if cap < n:
                    assert binomial_prefix(n)[cap+2] << j > 1 << power
    x = (1 << 1000)-1
    encoded = encode(x)
    assert (encoded >> M) & ((1 << 1000)-1) == x
    assert encoded.bit_length() <= 2000


def main():
    json.load(sys.stdin)
    fields = field_checks()
    table_checks()
    error, success, removal, occupied = certificate()
    manifest = json.loads((HERE/"submission.json").read_text())
    B,L = [
        Fraction(int(claim["score"]["numerator"]),int(claim["score"]["denominator"]))
        for claim in manifest["claims"]
    ]
    W,costs = operation_bound()
    assert error <= B and W**L.denominator <= 2**L.numerator
    print(json.dumps({
        "error_upper_bound":str(B),"ops":str(W),"log_ops_bound":str(L),
        "gcd_checks":fields,"cutoff":CUTOFF,"bands":BANDS,"memory_end":MEMORY_END,
        "success_before_single_removal_numerator":str(success),
        "weighted_removal_ceil_1e24":str(-(-removal.numerator*10**24//removal.denominator)),
        "last_occupied_lower_count":str(occupied[-1]),"operation_costs":costs,
    }))


if __name__ == "__main__":
    main()
