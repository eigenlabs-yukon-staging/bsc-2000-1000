"""Reproduce the exact coset spectra, then check the product bound."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    BLOCKS, COLUMNS, CORRECTION_BASE, LIST_BASE, LIST_END, LIST_SIZE, PARITY_BASE, PRIVATE_BITS, SIDE_BITS,
    TAG_BASE, TAG_BITS, WEIGHT_BASE, displacement, encode, operation_bound, quotient,
)

CUTOFF = 17
LIST_CUTOFF = 22
DIGIT_BITS = 51
DIGIT_MASK = (1 << DIGIT_BITS) - 1
SPECTRUM_MODULUS = 1 << ((2 * CUTOFF + 1) * DIGIT_BITS)


def gray_order(columns):
    values = [0]
    for column in columns:
        values.extend(value ^ column for value in reversed(values))
    return values


def pack(sequence, size):
    return sum(coefficient << (size * degree)
               for degree, coefficient in enumerate(sequence))


def coefficient(value, index, size):
    return (value >> (index * size)) & ((1 << size) - 1)


def side_spectrum(sector):
    # A word lists the private bits first, followed by the parity bits.
    basis = tuple((1 << i) | (col << PRIVATE_BITS) for i, col in enumerate(COLUMNS[:PRIVATE_BITS]))
    words = gray_order(basis)
    total = 0
    local = None
    parent = [0] * (SIDE_BITS + 1)
    for tag in range(sector << 8, (sector + 1) << 8):
        shift = displacement(tag) << PRIVATE_BITS
        hist = [0] * (SIDE_BITS + 1)
        for word in words:
            hist[(word ^ shift).bit_count()] += 1
        for i, count in enumerate(hist):
            parent[i] += count
        if tag == 0:
            local = hist
        else:
            packed = pack(hist[:2 * CUTOFF + 1], DIGIT_BITS)
            total = (total + packed * packed) % SPECTRUM_MODULUS
    return local, [coefficient(total, i, DIGIT_BITS) for i in range(2 * CUTOFF + 1)], parent


def cumulative_binomial(width):
    prefix = []
    for n in range(width + 1):
        row = []
        for w in range(n + 1):
            row.append((row[-1] if row else 0) + comb(n, w))
        prefix.append(row)
    return prefix


def pair_bounds(width, spectrum, cutoff=CUTOFF):
    prefix = cumulative_binomial(width)
    totals = []
    for t in range(cutoff + 1):
        pairs = 0
        for d in range(1, 2 * t + 1):
            common = 0
            for a in range(max(0, d - t), min(d, t) + 1):
                rest = min(width - d, t - max(a, d - a))
                common += comb(d, a) * prefix[width - d][rest]
            pairs += common * spectrum[d]
        totals.append(pairs)
    return totals


def coverage(volumes, pair_limits):
    cumulative = []
    previous = 0
    for volume, pairs in zip(volumes, pair_limits, strict=True):
        k = 1 + pairs // volume
        denom = k * (k + 1)
        lower = (2 * k * volume - pairs + denom - 1) // denom
        cauchy = (volume * volume + volume + pairs - 1) // (volume + pairs)
        previous = max(previous, lower, cauchy)
        cumulative.append(previous)
    return cumulative


def component_numerator(local_spectrum, outer_spectrum):
    local = coverage(cumulative_binomial(SIDE_BITS)[SIDE_BITS][:CUTOFF + 1],
                     pair_bounds(SIDE_BITS, local_spectrum))
    assert local[-1] <= 1 << 30
    digits = [v - (local[i - 1] if i else 0) for i, v in enumerate(local)]
    squared = pack(digits, 61) ** 2
    volumes = []
    for t in range(CUTOFF + 1):
        volumes.append((volumes[-1] if volumes else 0) + coefficient(squared, t, 61))
    occupied = coverage(volumes, pair_bounds(2 * SIDE_BITS, outer_spectrum))
    assert occupied[-1] <= 1 << 50
    numerator = sum((value - (occupied[t - 1] if t else 0)) * 9 ** (100 - t)
                    for t, value in enumerate(occupied))
    return numerator, occupied


def list_success_numerator(parent_spectrum):
    volumes = cumulative_binomial(SIDE_BITS)[SIDE_BITS]
    pair_limits = pair_bounds(SIDE_BITS, parent_spectrum, LIST_CUTOFF)
    previous = 0
    numerator = 0
    for t, pairs in enumerate(pair_limits):
        volume = volumes[t]
        k = max(2 * LIST_SIZE - 1, 1 + pairs // volume)
        denom = k * (k + 1)
        lower = (LIST_SIZE * (2 * k * volume - pairs) + denom - 1) // denom
        current = max(previous, lower)
        assert current <= LIST_SIZE * (1 << PRIVATE_BITS)
        numerator += (current - previous) * 9 ** (SIDE_BITS - t)
        previous = current
    return numerator


def product_error(local, outer, parent):
    numerator, occupied = component_numerator(local, outer)
    other = list_success_numerator(parent)
    full_success = Fraction(numerator, 10**100)
    local_list_success = Fraction(other, 10**50)
    retained = full_success - (1 - local_list_success)**2
    assert retained > 0
    return 1 - retained**BLOCKS, occupied, numerator, other


def main():
    sector = json.load(sys.stdin)["sector"]
    data = json.loads((HERE / "spectra.json").read_text())
    observed_local, observed_outer, observed_parent = side_spectrum(sector)
    assert observed_outer == data["sectors"][sector]
    assert observed_parent == data["parent_sectors"][sector]
    if observed_local is not None:
        assert observed_local == data["local"]
    for t in range(1 << TAG_BITS):
        assert displacement(t) & ((1 << TAG_BITS) - 1) == t
        assert quotient(displacement(t)) == 0
    for shift in range(1 << PRIVATE_BITS):
        assert quotient(shift << TAG_BITS) == shift
    assert all(0 <= column < 1 << 30 for column in COLUMNS)
    assert (0 < LIST_BASE < LIST_END <= WEIGHT_BASE < CORRECTION_BASE < PARITY_BASE < TAG_BASE
            < TAG_BASE + (1 << TAG_BITS) <= 1 << 32)
    assert encode((1 << 1000) - 1).bit_length() <= 2000
    spectrum = [sum(parts) for parts in zip(*data["sectors"], strict=True)]
    parent = [sum(parts) for parts in zip(*data["parent_sectors"], strict=True)]
    actual, occupied, numerator, list_numerator = product_error(data["local"], spectrum, parent)
    manifest = json.loads((HERE / "submission.json").read_text())
    error_bound, log_bound = [
        Fraction(int(claim["score"]["numerator"]), int(claim["score"]["denominator"]))
        for claim in manifest["claims"]
    ]
    ops = operation_bound()
    assert actual <= error_bound and ops ** log_bound.denominator <= 2 ** log_bound.numerator
    print(json.dumps({
        "sector": sector, "tag_count": 256, "private_count": 1 << PRIVATE_BITS,
        "checked_outer_spectrum": observed_outer,
        "checked_private_spectrum": observed_local,
        "checked_parent_sector": observed_parent,
        "occupied_syndromes_lower_bounds": occupied,
        "unrestricted_component_success_numerator": str(numerator),
        "side_list_success_numerator": str(list_numerator),
        "error_upper_bound": str(error_bound), "ops": ops, "log_ops_bound": str(log_bound),
    }))


if __name__ == "__main__":
    main()
