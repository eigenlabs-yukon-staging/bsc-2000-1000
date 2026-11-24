"""Exact finite checks for the three unequal component positions."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import CHUNK_BITS, F333, F334, FIELDS, RADIUS, encode, operation_bound, polynomial_remainder

BALANCE = 24
CUTOFF = 120


def gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def irreducibility(field):
    x = 2
    checks = {}
    for i in range(1, field.width + 1):
        x = field.square(x)
        if i in field.rabin_degrees:
            checks[i] = gcd(x ^ 2, field.modulus)
            assert checks[i] == 1
    assert x == 2
    return checks


def success_numerator(width):
    counts = [comb(width, i) for i in range(CUTOFF + 1)]
    balls = []
    bias = []
    for weight in range(CUTOFF + 1):
        balls.append(counts[weight] + (balls[-1] if balls else 0))
        coefficient = sum(comb(BALANCE, a) * comb(width - 1 - 2 * BALANCE, weight - 2 * a)
                          for a in range(min(BALANCE, weight // 2) + 1))
        assert coefficient <= balls[-1]
        bias.append(coefficient)
    bad = 3 * (balls[BALANCE] - 1 + balls[BALANCE - 1])
    size = 1 << width
    previous = total = 0
    for weight in range(CUTOFF + 1):
        volume = collision = 0
        for systematic in range(min(RADIUS, weight) + 1):
            count = counts[systematic]
            v, c = balls[weight - systematic], bias[weight - systematic]
            volume += count * v
            collision += count * (v + c) + min(count, bad) * (v - c)
        k = 1 + collision // size
        lower = max(previous, volume * size // (size + collision),
                    volume * (2 * k * size - collision) // (size * k * (k + 1)))
        assert previous <= lower <= size
        total += (lower - previous) * 9 ** (2 * width - weight)
        previous = lower
    return total


def test_tables(field):
    rng = random.Random(9120 + field.width)
    mask = (1 << field.width) - 1
    values = [0, 1, mask, 1 << (field.width - 1)] + [rng.getrandbits(field.width) for _ in range(6)]
    reduction = field.reduction_table(4)
    assert max(reduction).bit_length() <= 64 and field.high_bits >= 4
    assert field.chunks * CHUNK_BITS <= field.words * 64
    for value in values:
        squared = 0
        for j in range(field.bytes):
            squared ^= field.square_entry(j, (value >> (8 * j)) & 255)
        assert squared == field.square(value)
        times = [0, squared]
        for i in range(1, 4):
            new = field.remainder(squared << i)
            times.extend(old ^ new for old in times[:])
        result = 0
        for j in range(field.nibbles - 1, -1, -1):
            result = ((result << 4) & mask) ^ reduction[result >> (field.width - 4)]
            result ^= times[(value >> (4 * j)) & 15]
        assert result == field.cube(value)
    indices = (0, 14, 15, 63, 64, 127, 128, field.width - 1)
    for direction in indices:
        value = rng.getrandbits(field.width)
        delta = 0
        for chunk in range(field.chunks):
            key = (value >> (CHUNK_BITS * chunk)) & ((1 << CHUNK_BITS) - 1)
            delta ^= field.affine_entry(direction, chunk, key)
        assert field.cube(value) ^ delta == field.cube(value ^ (1 << direction))
    return len(values), len(indices)


def main():
    json.load(sys.stdin)
    md = json.loads((HERE / "submission.json").read_text())
    error_bound, log_bound = [
        Fraction(int(claim["score"]["numerator"]), int(claim["score"]["denominator"]))
        for claim in md["claims"]
    ]
    checks = {}
    for field in (F333, F334):
        checks[field.width] = {
            "gcd_checks": irreducibility(field),
            "success_numerator": str(success_numerator(field.width)),
            "bank": [field.start, field.end],
            "arithmetic_samples": test_tables(field),
            "cost_data": field.operation_bound()[1],
        }
    correct = int(checks[333]["success_numerator"]) ** 2 * int(checks[334]["success_numerator"])
    actual = 1 - Fraction(correct, 10 ** 2000)
    W = operation_bound()
    assert F333.start == 48 <= F333.end == F334.start <= F334.end <= 1 << 32
    assert sum(field.width for field in FIELDS) == 1000
    word = encode((1 << 1000) - 1)
    for position, field in zip((0, 666, 1332), FIELDS):
        assert (word >> position) & ((1 << field.width) - 1) == (1 << field.width) - 1
    assert word.bit_length() <= 2000
    assert actual <= error_bound and W ** log_bound.denominator <= 2 ** log_bound.numerator
    print(json.dumps({
        "fields": checks, "error_upper_bound": str(error_bound),
        "ops_upper_bound": str(W), "log_ops_bound": str(log_bound),
    }))


if __name__ == "__main__":
    main()
