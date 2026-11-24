"""Exact finite checks for the uniform component layout."""
from bisect import bisect_left
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import CHUNK_BITS, EVEN, FIELD, FIELDS, MAPS, OMEGA, PARITY_RADIUS, RADIUS, ROOT_PROGRAM, carryless, encode, operation_bound, polynomial_remainder

CUTOFF = 90


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


def imbalance_envelopes(width, counts):
    rows = [[] for _ in counts]
    n = width - 1
    for weight in range(1, width + 1):
        multiplicity = 3 * comb(width, weight)
        gap = n - 2 * (weight - 1)
        previous, signed = 0, 1
        for radius, row in enumerate(rows):
            row.append((abs(signed), multiplicity))
            if radius < CUTOFF:
                successor, remainder = divmod(
                    gap * signed - (n - radius + 1) * previous, radius + 1
                )
                assert remainder == 0
                previous, signed = signed, successor
    envelopes = []
    for row in rows:
        ordered = sorted(row, reverse=True)
        caps = [0]
        totals = [0]
        for value, count in ordered:
            caps.append(caps[-1] + count)
            totals.append(totals[-1] + value * count)
        envelope = []
        for capacity in counts:
            i = bisect_left(caps, capacity) - 1
            envelope.append(totals[i] + (capacity - caps[i]) * ordered[i][0])
        envelopes.append(envelope)
    return envelopes


def success_numerator(width):
    counts = [comb(width, i) for i in range(CUTOFF + 1)]
    balls = []
    for count in counts:
        balls.append(count + (balls[-1] if balls else 0))
    imbalance = imbalance_envelopes(width, counts)
    size = 1 << width
    previous = total = 0
    for weight in range(CUTOFF + 1):
        volume = collision = 0
        for systematic in range(weight + 1):
            radius = weight - systematic
            if systematic > RADIUS:
                radius = min(radius, PARITY_RADIUS)
            term = counts[systematic] * balls[radius]
            volume += term
            collision += term + imbalance[radius][systematic]
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
    indices = (0, CHUNK_BITS - 1, CHUNK_BITS, 63, 64, field.width - 1)
    for direction in indices:
        value = rng.getrandbits(field.width)
        delta = 0
        for chunk in range(field.chunks):
            key = (value >> (CHUNK_BITS * chunk)) & ((1 << CHUNK_BITS) - 1)
            delta ^= field.affine_entry(direction, chunk, key)
        assert field.cube(value) ^ delta == field.cube(value ^ (1 << direction))
    return len(values), len(indices)


def root_circuit(field, value):
    arrays = {"KEY": value}
    for instruction in ROOT_PROGRAM:
        op, *parameters = instruction
        if op == "COPY":
            source, dest = parameters
            arrays[dest] = arrays[source]
        elif op == "MUL":
            a, b, dest = parameters
            arrays[dest] = field.remainder(carryless(arrays[a], arrays[b]))
        else:
            tag, source, dest = parameters
            x = arrays[source]
            result = 0
            for j in range(field.map_slices):
                key = (x >> (16 * j)) & 65535
                result ^= field.map_entry(tag, j, key)
            arrays[dest] = result
    return arrays["Z"]


def root_checks():
    # The same straight-line macros, interpreted on exponents instead of words.
    exponents = {"KEY": 1}
    for instruction in ROOT_PROGRAM:
        op, *parameters = instruction
        if op == "COPY":
            source, dest = parameters
            exponents[dest] = exponents[source]
        elif op == "MUL":
            a, b, dest = parameters
            exponents[dest] = exponents[a] + exponents[b]
        else:
            tag, source, dest = parameters
            exponents[dest] = exponents[source] << tag
    q = 1 << FIELD.width
    e = exponents["Z"]
    if EVEN:
        third = (q-1)//3
        assert 3*e % third == 1 and FIELD.square(OMEGA) == (OMEGA ^ 1)
        assert OMEGA not in (0,1)
    else:
        assert 3*e % (q-1) == 1
    rng = random.Random(882511)
    values = [0, 1] + [rng.getrandbits(FIELD.width) for _ in range(6)]
    for value in values:
        x = FIELD.cube(value)
        a = root_circuit(FIELD, x)
        assert a == FIELD.power(x, e) and FIELD.cube(a) == x
        if EVEN:
            b = 0
            for j in range(FIELD.map_slices):
                b ^= FIELD.map_entry("omega",j,(a>>(16*j))&65535)
            assert FIELD.cube(b) == FIELD.cube(a ^ b) == x
            assert value in (a, b, a ^ b)
        else:
            assert a == value
    return len(values), str(e), len(MAPS)


def main():
    json.load(sys.stdin)
    md = json.loads((HERE / "submission.json").read_text())
    error_bound, log_bound = [
        Fraction(int(claim["score"]["numerator"]), int(claim["score"]["denominator"]))
        for claim in md["claims"]
    ]
    gcd_checks = irreducibility(FIELD)
    numerator = success_numerator(FIELD.width)
    actual = 1 - Fraction(numerator ** len(FIELDS), 10 ** 2000)
    _, cost_data = FIELD.operation_bound()
    W = operation_bound()
    assert FIELD.start == 48 < FIELD.end <= 1 << 32
    assert sum(field.width for field in FIELDS) == 1000
    word = encode((1 << 1000) - 1)
    for i,field in enumerate(FIELDS):
        assert (word >> (2 * FIELD.width * i)) & ((1 << field.width) - 1) == (1 << field.width) - 1
    assert word.bit_length() <= 2000
    assert actual <= error_bound and W ** log_bound.denominator <= 2 ** log_bound.numerator
    print(json.dumps({
        "gcd_checks": gcd_checks, "component_success_numerator": str(numerator),
        "envelope_radii": CUTOFF+1, "mask_weights": FIELD.width,
        "memory_end": FIELD.end, "cost_data": cost_data,
        "arithmetic_samples": test_tables(FIELD), "root_arithmetic": root_checks(),
        "error_upper_bound": str(error_bound),
        "ops_upper_bound": str(W), "log_ops_bound": str(log_bound),
    }))


if __name__ == "__main__":
    main()
