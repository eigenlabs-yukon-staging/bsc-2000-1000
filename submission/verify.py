"""Exact field and received-word occupancy checks for the fixed product code."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    COMPONENTS, HIGH_BITS, MEMORY_END, MODULUS, NIBBLES, RADIUS, WIDTH, cube,
    operation_bound, reduction_table, remainder, square,
)

BALANCE = 24
CUTOFF = 165


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def irreducibility():
    x = 2
    checks = {}
    for i in range(1, WIDTH + 1):
        x = square(x)
        if i in (WIDTH // 2, WIDTH // 5):
            checks[i] = gcd(x ^ 2, MODULUS)
            assert checks[i] == 1
    assert x == 2
    return checks


def success_numerator():
    m = WIDTH
    counts = [comb(m, i) for i in range(CUTOFF + 1)]
    balls = []
    bias = []
    for t in range(CUTOFF + 1):
        balls.append(counts[t] + (balls[-1] if balls else 0))
        coefficient = sum(comb(BALANCE, a) * comb(m - 1 - 2 * BALANCE, t - 2 * a)
                          for a in range(min(BALANCE, t // 2) + 1))
        assert coefficient <= balls[-1]
        bias.append(coefficient)
    bad = 3 * (balls[BALANCE] - 1 + balls[BALANCE - 1])
    Q = 1 << m
    prev = total = 0
    for t in range(CUTOFF + 1):
        volume = collision = 0
        for s in range(min(RADIUS, t) + 1):
            n = counts[s]
            v, c = balls[t - s], bias[t - s]
            volume += n * v
            collision += n * (v + c) + min(n, bad) * (v - c)
        k = 1 + collision // Q
        lower = max(prev, volume * Q // (Q + collision),
                    volume * (2 * k * Q - collision) // (Q * k * (k + 1)))
        assert prev <= lower <= Q
        total += (lower - prev) * 9 ** (2 * m - t)
        prev = lower
    return total


def test_window_arithmetic():
    rng = random.Random(9120)
    mask = (1 << WIDTH) - 1
    inputs = [0, 1, mask, 1 << (WIDTH - 1)] + [rng.getrandbits(WIDTH) for _ in range(8)]
    reduction = reduction_table(4)
    assert max(reduction).bit_length() <= 64 and HIGH_BITS >= 4
    for a in inputs:
        operands = [0, square(a)]
        for i in range(1, 4):
            new = remainder(operands[1] << i)
            operands.extend(x ^ new for x in operands[:])
        result = 0
        for j in range(NIBBLES - 1, -1, -1):
            result = ((result << 4) & mask) ^ reduction[result >> (WIDTH - 4)]
            result ^= operands[(a >> (4 * j)) & 15]
        assert result == cube(a)
    return len(inputs)


def main():
    json.load(sys.stdin)
    md = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in md["claims"]]
    checks = irreducibility()
    success = success_numerator()
    actual = 1 - Fraction(success ** COMPONENTS, 10 ** 2000)
    nodes, W, macro_costs = operation_bound()
    assert WIDTH * COMPONENTS == 1000 and MEMORY_END <= 1 << 32
    assert actual <= B and W ** L.denominator <= 2 ** L.numerator
    print(json.dumps({
        "gcd_checks": checks, "component_success_numerator": str(success),
        "error_upper_bound": str(B), "ops_upper_bound": str(W),
        "log_ops_bound": str(L), "candidates_per_component": str(nodes),
        "square_mul_score_counts": macro_costs, "memory_end_exclusive": MEMORY_END,
        "window_sample_count": test_window_arithmetic(),
    }))


if __name__ == "__main__":
    main()
