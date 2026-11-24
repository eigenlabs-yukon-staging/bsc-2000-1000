from math import comb
from fractions import Fraction
from collections import Counter
from pathlib import Path
import json
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    AS_map, AS_pivots, CUBE_EXP, LAMBDA, MASK, MEMORY_END, MODULUS, PAIRS, RADIUS,
    THIRD, W_bound, chart, codeword, constant, cube, cube_root, inverse, map_tags,
    power, product, reduction_table, remainder, square,
)

WIDTH=500
CUTOFF=240
HALF=1<<WIDTH
Q=HALF**2
GOOD=(HALF-1)**2
COEFFICIENTS=[comb(WIDTH,w) for w in range(CUTOFF+1)]

def prefix(items):
    ans=[0]
    for x in items:ans.append(ans[-1]+x)
    return ans

def overlapping_tail_sum(pref,weight,radius):
    if radius<0:return 0
    if 2*radius+1>=weight:return pref[-1]
    return pref[radius+1]+pref[-1]-pref[weight-radius]

def parity_profiles(r,D):
    c=COEFFICIENTS
    bounded=[sum(comb(D,j)*comb(WIDTH-2*D,w-2*j) for j in range(min(D,w//2)+1))
             for w in range(CUTOFF+1)]
    shell_vol=[];shell_bias=[]
    for w in range(CUTOFF+1):
        shell_vol.append(prefix(c[p]*c[w-p] for p in range(w+1)))
        shell_bias.append(prefix(bounded[p]*c[w-p] for p in range(w+1)))
    candidates={}
    for a in range(-1,r+1):
        vals=[]
        v=u=0
        for w in range(CUTOFF+1):
            if w<=r:
                v+=shell_vol[w][-1];u+=shell_bias[w][-1]
            else:
                v+=overlapping_tail_sum(shell_vol[w],w,a)
                u+=overlapping_tail_sum(shell_bias[w],w,a)
            vals.append((v,u))
        candidates[a]=vals
    v=u=0; full=[]
    for w in range(CUTOFF+1):
        v+=shell_vol[w][-1];u+=shell_bias[w][-1]
        full.append((v,u))
    return candidates,full


def certificate(r,D):
    parity,complete=parity_profiles(r,D)
    bad=15*(2*sum(COEFFICIENTS[:D]))**2
    sys=[]
    for s in range(CUTOFF+1):
        if s<=r:
            sys.append([(complete,comb(2*WIDTH,s))])
        else:
            totals=Counter()
            for a in range(s//2+1):
                count=COEFFICIENTS[a]*COEFFICIENTS[s-a]
                if 2*a<s: count*=2
                totals[max(-1,r-a)]+=count
            sys.append([(parity[a], count) for a,count in totals.items()])
    success=previous=0
    for t in range(CUTOFF+1):
        V=A=0
        for s in range(t+1):
            remain=t-s
            for shape,n in sys[s]:
                v,b=shape[remain]
                V+=n*v
                A+=n*(v+b)+min(n,bad)*(v-b)
        d=1+A//GOOD
        low=max(previous,V*GOOD**2//(Q*(GOOD+A)),V*(2*d*GOOD-A)//(Q*d*(d+1)))
        assert previous<=low<=Q
        success+=(low-previous)*9**(2000-t)
        previous=low
    return 1-Fraction(success,10**2000)

def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def check_field():
    x = 2
    out = {}
    for i in range(1, WIDTH + 1):
        x = square(x)
        if i in (100, 250):
            out[i] = gcd(x ^ 2, MODULUS)
            assert out[i] == 1
    assert x == 2
    assert 3 * THIRD == MASK and 3 * CUBE_EXP == 2 * THIRD + 1
    assert CUBE_EXP == 1 + 56 * sum(1 << (6 * i) for i in range(83))
    zeta = constant("zeta")
    assert zeta != 1 and cube(zeta) == 1 and zeta == power(LAMBDA, THIRD)
    assert product(LAMBDA, constant("lambda_inv")) == 1
    assert max(reduction_table(4)).bit_length() <= 64
    pivots = AS_pivots()
    assert len(pivots) == WIDTH - 1
    for _, image, _ in pivots:
        x = AS_map(image)
        assert square(x) ^ x == image
    return out


def check_charts():
    rng = random.Random(266310)
    inputs = [(1, 1), (MASK, 1)]
    inputs += [(rng.randrange(1, HALF), rng.randrange(1, HALF)) for _ in range(3)]
    for message in inputs:
        a, b = message
        word = codeword(a, b)
        assert product(a, inverse(a)) == 1 and cube(cube_root(cube(b))) == cube(b)
        for i, j in PAIRS:
            candidates = chart((i, j), word[i], word[j])
            assert (a, b) in candidates and len(candidates) <= 6
            assert all(a0 and b0 and (codeword(a0, b0)[i], codeword(a0, b0)[j]) ==
                       (word[i], word[j]) for a0, b0 in candidates)
    assert not chart((2, 3), 0, 1)
    assert not chart((0, 2), 0, 7)
    return len(inputs) * len(PAIRS)


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in manifest["claims"]]
    field = check_field()
    n, W = W_bound()
    b = certificate(RADIUS, 42)
    assert b <= B and W ** L.denominator <= 2 ** L.numerator
    assert MEMORY_END <= 1 << 32
    print(json.dumps({"gcd_checks": field, "third_root_hex": hex(constant("zeta")),
                      "exact_error_bound": str(B), "error_at_21_decimal_places_ceiling":
                          -(-(b.numerator * 10 ** 21) // b.denominator),
                      "map_count": len(map_tags()), "memory_end_exclusive": MEMORY_END,
                      "search_nodes_per_chart": str(n), "W": str(W), "L": str(L),
                      "chart_sample_count": check_charts()}))


if __name__ == "__main__":
    main()
