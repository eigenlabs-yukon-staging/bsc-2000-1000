"""Exact field checks and finite-length occupancy certificate."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    CHAIN, EXPONENT, MASK, MEMORY_END, MODULUS, RADIUS_FORWARD, RADIUS_REVERSE, THIRD, WIDTH,
    cube, inverse_chain, map_tags, operation_bound, power, product,
    reduction_table, remainder, root_of_unity, square,
)

BALANCE = 48
CUTOFF = 240


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


def profiles():
    m = WIDTH
    combos = [comb(m, t) for t in range(CUTOFF + 1)]
    balls = []
    bias = []
    for t in range(CUTOFF + 1):
        balls.append(combos[t] + (balls[-1] if balls else 0))
        c = sum(comb(BALANCE, k) * comb(m - 1 - 2 * BALANCE, t - 2 * k)
                for k in range(min(BALANCE, t // 2) + 1))
        assert c <= balls[t]
        bias.append(c)
    nbad = 3 * (balls[BALANCE] - 1 + balls[BALANCE - 1])
    for t in range(CUTOFF + 1):
        volume = collisions_per_volume_num = 0
        for d in range(t + 1):
            ell = t - d
            if d > RADIUS_FORWARD:
                ell = min(ell, RADIUS_REVERSE)
            count = combos[d]
            v, c = balls[ell], bias[ell]
            volume += count * v
            collisions_per_volume_num += count * (v + c) + min(count, nbad) * (v - c)
        yield volume, collisions_per_volume_num


def error_bound():
    scale = 1 << WIDTH
    prev = numerator = 0
    checkpoints = {}
    for t, (volume, collisions_per_volume_num) in enumerate(profiles()):
        a = collisions_per_volume_num
        d = 1 + a // scale
        lower = max(prev, volume * scale // (scale + a),
                    volume * (2 * d * scale - a) // (scale * d * (d + 1)))
        assert prev <= lower <= scale
        numerator += (lower - prev) * 9 ** (2 * WIDTH - t)
        prev = lower
        if 220 <= t <= 225:
            checkpoints[t] = str(lower)
    return 1 - Fraction(numerator, 10 ** (2 * WIDTH)), checkpoints


def check_chain():
    assert THIRD * 3 == (1 << WIDTH) - 1 and 3 * EXPONENT == THIRD + 1
    length = 1
    for bit in CHAIN:
        length = 2 * length + bit
    assert length == 166
    assert EXPONENT == 2 + 112 * sum(1 << (6 * i) for i in range(length))
    zeta = root_of_unity()
    assert 0 < zeta <= MASK and zeta != 1 and cube(zeta) == 1
    rng = random.Random(490325)
    inputs = [0, 1, 7, 1 << (WIDTH - 1)] + [rng.getrandbits(WIDTH) for _ in range(6)]
    for z in inputs:
        root = inverse_chain(z)
        assert root == power(z, EXPONENT)
        if cube(root) == z:
            assert cube(product(root, zeta)) == z
            assert cube(product(root, square(zeta))) == z
        assert cube(inverse_chain(cube(z))) == cube(z)
    assert max(reduction_table(4)).bit_length() <= 64
    return len(inputs), hex(zeta)


def main():
    json.load(sys.stdin)
    md = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in md["claims"]]
    checks = verify_modulus()
    samples, zeta = check_chain()
    bound, checkpoints = error_bound()
    candidates, W = operation_bound()
    assert bound <= B and W ** L.denominator <= 2 ** L.numerator
    assert MEMORY_END <= 1 << 32
    rounded = -(-(bound.numerator * 10 ** 21) // bound.denominator)
    print(json.dumps({
        "irreducibility_gcd_checks": checks, "third_root_of_unity": zeta,
        "chain_sample_count": samples, "error_upper_bound": str(B),
        "raw_bound_ceil_at_21_decimal_places": rounded,
        "coverage_checkpoints": checkpoints, "candidate_counts": [str(n) for n in candidates],
        "W": str(W), "L": str(L), "map_tags": map_tags(), "memory_end_exclusive": MEMORY_END,
    }))


if __name__ == "__main__":
    main()
