"""Exact irreducibility, rational error, and worst-case exponent checks."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    MASK, MEMORY_END, MODULUS, RADIUS, WIDTH, cube, operation_bound,
    reduction_table, remainder, square,
)

BALANCE = 48
NOISE_CUTOFF = 222


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def verify_modulus():
    x = 2
    checks = {}
    for i in range(1, WIDTH + 1):
        x = square(x)
        if i in (WIDTH // 5, WIDTH // 2):
            checks[i] = gcd(x ^ 2, MODULUS)
            assert checks[i] == 1
    assert x == 2
    return checks


def cumulative(values):
    total = 0
    for v in values:
        total += v
        yield total


def error_bound():
    m = WIDTH
    combinations = [comb(m, i) for i in range(m + 1)]
    balls = list(cumulative(combinations))
    masses = [c * 9 ** (m - i) for i, c in enumerate(combinations)]
    tails = [0] * (m + 2)
    for i in range(m, -1, -1):
        tails[i] = tails[i + 1] + masses[i]
    assert tails[0] == 10 ** m
    bias = []
    for ell in range(NOISE_CUTOFF):
        bound = sum(comb(BALANCE, k) * comb(m - 1 - 2 * BALANCE, ell - 2 * k)
                    for k in range(min(BALANCE, ell // 2) + 1))
        assert 0 <= bound <= balls[ell]
        bias.append(bound)
    bad = 3 * (balls[BALANCE] - 1 + balls[BALANCE - 1])
    scale = 1 << (m + 1)
    conditional = []
    for w in range(NOISE_CUTOFF):
        total = 0
        for t in range(min(w, RADIUS) + 1):
            ell = w - t
            volume = balls[ell] + (balls[ell - 1] if ell else 0)
            sign = bias[ell] + (bias[ell - 1] if ell else 0)
            count = combinations[t]
            total += count * (volume + sign) + min(count, bad) * (volume - sign)
        conditional.append(min(total, scale))
    numerator = tails[RADIUS + 1] * 10 ** m * scale
    for i in range(RADIUS + 1):
        cutoff = max(0, NOISE_CUTOFF - i)
        subtotal = tails[cutoff] * scale
        subtotal += sum(masses[j] * conditional[i + j] for j in range(cutoff))
        numerator += masses[i] * subtotal
    return Fraction(numerator, scale * 10 ** (2 * m)), bad


def mul_model(a, b):
    """Host illustration of the RAM window identity, word splitting omitted."""
    red = reduction_table(4)
    table = [0, a]
    for j in range(1, 4):
        direction = remainder(a << j)
        table += [v ^ direction for v in table]
    r = 0
    for j in range(249, -1, -1):
        r = ((r << 4) & MASK) ^ red[r >> (WIDTH - 4)] ^ table[(b >> (4 * j)) & 15]
    return r


def main():
    json.load(sys.stdin)
    md = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in md["claims"]]
    checks = verify_modulus()
    actual_bound, bad = error_bound()
    candidates, W = operation_bound()
    assert actual_bound <= B
    assert W ** L.denominator <= 2 ** L.numerator
    assert MEMORY_END <= 1 << 32
    rng = random.Random(136420)
    inputs = [0, 1, MASK, 1 << (WIDTH - 1)] + [rng.getrandbits(WIDTH) for _ in range(8)]
    assert all(mul_model(square(a), a) == cube(a) for a in inputs)
    rounded = -(-(actual_bound.numerator * (10 ** 21)) // actual_bound.denominator)
    print(json.dumps({
        "irreducibility_gcd_checks": checks, "final_frobenius_is_x": True,
        "error_upper_bound": str(B), "raw_bound_ceil_at_21_decimal_places": rounded,
        "max_unbalanced_differences": str(bad), "candidate_count": str(candidates),
        "W": str(W), "L": str(L), "memory_end_exclusive": MEMORY_END,
        "window_identity_sample_count": len(inputs),
    }))


if __name__ == "__main__":
    main()
