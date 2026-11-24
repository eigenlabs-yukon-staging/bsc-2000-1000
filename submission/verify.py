from math import comb
from functools import lru_cache
from fractions import Fraction
from collections import Counter
from pathlib import Path
import json
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    AS_map, AS_pivots, CUBE_EXP, LAMBDA, MASK, MEMORY_END, MODULUS, PAIRS, RADII,
    THIRD, W_bound, chart, codeword, constant, cube, cube_root, derivative_words, inverse, map_tags,
    power, product, reduction_table, remainder, square,
)

M=500
WIDTH=M
HALF=1<<M
CUT=240
Q=1<<1000
G=((1<<500)-1)**2
D=42
C=[comb(M,w) for w in range(CUT+1)]
CC=[comb(2*M,w) for w in range(CUT+1)]
CBOUNDS=[sum(comb(D,j)*comb(M-2*D,i-2*j) for j in range(min(D,i//2)+1))
         for i in range(CUT+1)]
BAD=15*(2*sum(C[:D]))**2


def prefix(seq):
    out=[0]
    for z in seq:out.append(out[-1]+z)
    return out


def region(pref,w,a,b):
    if a+b+1>=w:
        return pref[-1]
    return pref[a+1]+pref[-1]-pref[w-b]


@lru_cache(1)
def shells():
    va=[]
    a1=[]
    a2=[]
    for w in range(CUT+1):
        va.append(prefix(C[j]*C[w-j] for j in range(w+1)))
        a1.append(prefix(CBOUNDS[j]*C[w-j] for j in range(w+1)))
        a2.append(prefix(C[j]*CBOUNDS[w-j] for j in range(w+1)))
    full=[]
    v=c=0
    for w in range(CUT+1):
        v+=va[w][-1]
        c+=a1[w][-1]
        full.append((v,c))
    return va,a1,a2,full


@lru_cache(4)
def profiles(rmul,rcube,rpar):
    va,ax,bx,full=shells()
    shapes={}
    for smallest in range(max(rmul,rcube)+2):
        a=max(-1,rmul-smallest)
        b=max(-1,rcube-smallest)
        key=(a,b)
        if key in shapes:continue
        v=x=y=0
        p=[]
        for w in range(CUT+1):
            if w<=rpar:
                v+=va[w][-1]
                x+=ax[w][-1]
                y+=bx[w][-1]
            else:
                v+=region(va[w],w,a,b)
                x+=region(ax[w],w,a,b)
                y+=region(bx[w],w,a,b)
            p.append((v,max(x,y)))
        shapes[key]=p
    return full,shapes


def certificate(rsys,rmul,rcube,rpar):
    full,shapes=profiles(rmul,rcube,rpar)
    systematic=[]
    for total in range(CUT+1):
        if total<=rsys:
            systematic.append([(CC[total],full)])
        else:
            bins=Counter()
            for a in range(total//2+1):
                bins[(max(-1,rmul-a),max(-1,rcube-a))]+=C[a]*C[total-a]*(1 if 2*a==total else 2)
            systematic.append([(n,shapes[key]) for key,n in bins.items()])
    prior=success=0
    for t in range(CUT+1):
        V=A=0
        for i in range(t+1):
            left=t-i
            for count,shape in systematic[i]:
                v,b=shape[left]
                V+=count*v
                A+=count*(v+b)+min(count,BAD)*(v-b)
        k=1+A//G
        low=max(prior,V*G*G//(Q*(G+A)),V*(2*k*G-A)//(Q*k*(k+1)))
        assert prior <= low <= Q
        success+=(low-prior)*9**(2000-t)
        prior=low
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
            assert (a, b) in candidates and len(candidates) <= 3
            assert all(a0 and b0 and (codeword(a0, b0)[i], codeword(a0, b0)[j]) ==
                       (word[i], word[j]) for a0, b0 in candidates)
    assert not chart((2, 3), 1, 0)
    assert not chart((2, 3), 0, 1)
    assert not chart((0, 2), 0, 7)
    return len(inputs) * len(PAIRS)


def check_updates():
    rng = random.Random(9014)
    for sample, i in enumerate((0, 1, 63, 64, 127, 498, 499)):
        a, b = (rng.randrange(1, HALF) for _ in range(2))
        original = codeword(a, b)
        chart_type = sample % 3
        if chart_type == 0:
            u, alpha, scalar = b, LAMBDA, a
            expected = codeword(a, b ^ (1 << i))
            recovered_index = 2
        else:
            outer = a if chart_type == 1 else b
            scalar = inverse(outer)
            alpha = cube(scalar)
            if chart_type == 1:
                alpha = product(LAMBDA, alpha)
            u = original[2]
            c_new = u ^ (1 << i)
            pair = (a, product(c_new, scalar)) if chart_type == 1 else (product(c_new, scalar), b)
            expected = codeword(*pair)
            recovered_index = chart_type % 2
        d = 0
        desired = [(j * 256 + ((u >> (8 * j)) & 255)) * 8 for j in range(63)]
        words = iter(derivative_words(i, alpha))
        cursor = 0
        for offset in desired:
            for _ in range(offset - cursor):
                next(words)
            for k in range(8):
                d ^= next(words) << (64 * k)
            cursor = offset + 8
        assert original[3] ^ d == expected[3]
        direction = scalar
        for _ in range(i):
            direction = remainder(direction << 1)
        assert original[recovered_index] ^ direction == expected[recovered_index]
    return 7


def check_table_recurrence():
    alpha = 0x2651a3bb73d2
    f1 = f2 = f3 = alpha
    checked = 0
    for i in range(M):
        if i in (0, 63, 499):
            lin, quad = f1, f2
            expected = iter(derivative_words(i, alpha))
            for byte in range(63):
                entries = [f3 if byte == 0 else 0]
                for bit in range(8):
                    column = 0 if byte == 62 and bit >= 4 else lin ^ quad
                    entries += [row ^ column for row in entries]
                    lin = remainder(lin << 1)
                    quad = remainder(quad << 2)
                for entry in entries:
                    for word in range(8):
                        assert next(expected) == (entry >> (64 * word)) & ((1 << 64) - 1)
                checked += len(entries)
        f1 = remainder(f1 << 2)
        f2 = remainder(f2 << 1)
        f3 = remainder(f3 << 3)
    return checked


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in manifest["claims"]]
    field = check_field()
    n, o, W = W_bound()
    b = certificate(RADII[0], RADII[1], RADII[2], RADII[5])
    assert b <= B and W ** L.denominator <= 2 ** L.numerator
    assert MEMORY_END <= 1 << 32
    print(json.dumps({"gcd_checks": field, "third_root_hex": hex(constant("zeta")),
                      "exact_error_bound": str(B), "error_at_21_decimal_places_ceiling":
                          -(-(b.numerator * 10 ** 21) // b.denominator),
                      "map_count": len(map_tags()), "memory_end_exclusive": MEMORY_END,
                      "search_nodes_per_chart": [str(x) for x in n],
                      "outer_nodes_per_chart": [str(x) for x in o], "W": str(W), "L": str(L),
                      "chart_sample_count": check_charts(),
                      "update_sample_count": check_updates(),
                      "table_entries_checked": check_table_recurrence()}))


if __name__ == "__main__":
    main()
