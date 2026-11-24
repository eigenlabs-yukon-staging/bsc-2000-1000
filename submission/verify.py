"""Reproduce the fixed local spectra and the exact line-code certificate."""
import base64
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import random
import sys
import zlib

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    BLOCKS, COLUMNS, END, GAMMA_BASE, GROUPS, HIGH, PACKED_BASE, PARITY_BASE,
    PRIVATE_BITS, PROFILE_BASE, REPLICATE, SIDE_BITS, SLOTS, TAG_BASE, LEADER_BASE, directions,
    displacement, encode, gamma, operation_bound, packed_address, profile_address, tag_tuple,
)

CUTOFF = 25
SECTOR_BOUNDARIES = (0, 1248, 2496, 3744, 4096)
DIGIT_BITS = 78
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


def read_histograms():
    streams = []
    for i in range(4):
        data = json.loads((HERE / f"spectra{i}.json").read_text())
        raw_block = zlib.decompress(base64.b64decode(data["packed_spectra"], validate=True))
        assert len(raw_block) == (SLOTS // 4) * (SIDE_BITS + 1) * 3
        streams.append(raw_block)
    raw = b"".join(streams)
    rows = []
    for tag in range(SLOTS):
        start = (SIDE_BITS + 1) * 3 * tag
        rows.append([int.from_bytes(raw[start+3*i:start+3*i+3], "little") for i in range(SIDE_BITS+1)])
        assert sum(rows[-1]) == 1 << PRIVATE_BITS
    return rows

def check_sector(sector, hists):
    # A word lists the private bits first, followed by the parity bits.
    basis = tuple((1 << i) | (col << PRIVATE_BITS) for i, col in enumerate(COLUMNS[:PRIVATE_BITS]))
    words = gray_order(basis)
    for tag in range(SECTOR_BOUNDARIES[sector], SECTOR_BOUNDARIES[sector + 1]):
        shift = displacement(tag) << PRIVATE_BITS
        hist = [0] * (SIDE_BITS + 1)
        for word in words:
            hist[(word ^ shift).bit_count()] += 1
        assert hist == hists[tag]


def outer_spectra(hists):
    # Bound every retained coefficient even for the sum of all terms.
    evaluation_scale = 1 << SIDE_BITS
    max_evaluation = max(evaluation_scale, max(
        sum(value << (SIDE_BITS-i) for i,value in enumerate(row)) for row in hists
    ))
    assert SLOTS**2 * (1 << (2*CUTOFF)) * max_evaluation**4 < (1 << DIGIT_BITS) * evaluation_scale**4
    buckets = [pack(row, DIGIT_BITS) for row in hists]
    mask = SPECTRUM_MODULUS - 1
    four = sum(pow(value, 4, SPECTRUM_MODULUS) for value in buckets[1:]) & mask
    three = 0
    lines = 0
    subspaces = 0
    for s, d in directions():
        subspaces += 1
        three = (three + 12 * buckets[s] * buckets[d] * buckets[s ^ d]) & mask
        for a in range(1, SLOTS):
            if a > min(a ^ s, a ^ d, a ^ s ^ d):
                continue
            lines += 1
            x = (buckets[a] * buckets[a ^ s]) & mask
            y = (buckets[a ^ d] * buckets[a ^ s ^ d]) & mask
            four = (four + 12 * x * y) & mask
    assert subspaces == (SLOTS - 1) // 3
    assert lines == subspaces * (SLOTS // 4 - 1)
    assert 1 + (SLOTS - 1) + 12 * (lines + subspaces) == SLOTS ** 2
    return [[coefficient(poly, i, DIGIT_BITS) for i in range(2*CUTOFF+1)] for poly in (four, three)]


def cumulative_binomial(width):
    prefix = []
    for n in range(width + 1):
        row = []
        for w in range(n + 1):
            row.append((row[-1] if row else 0) + comb(n, w))
        prefix.append(row)
    return prefix


def pair_bounds(width, spectrum):
    prefix = cumulative_binomial(width)
    totals = []
    for t in range(CUTOFF + 1):
        pairs = 0
        for d in range(1, 2 * t + 1):
            if spectrum[d]:
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


def component_numerator(hists):
    local = coverage(cumulative_binomial(SIDE_BITS)[SIDE_BITS][:CUTOFF + 1],
                     pair_bounds(SIDE_BITS, hists[0]))
    assert local[-1] <= 1 << 31
    digits = [v - (local[i - 1] if i else 0) for i, v in enumerate(local)]
    fourth_power = pack(digits, 125) ** 4
    volumes = []
    for t in range(CUTOFF + 1):
        volumes.append((volumes[-1] if volumes else 0) + coefficient(fourth_power, t, 125))
    four, three = outer_spectra(hists)
    quad_pairs = pair_bounds(4 * SIDE_BITS, four)
    triple_pairs = pair_bounds(3 * SIDE_BITS, three)
    pairs = [
        quad_pairs[t] + sum(min(comb(SIDE_BITS, w), 1 << 31) * triple_pairs[t-w]
                            for w in range(t+1))
        for t in range(CUTOFF + 1)
    ]
    occupied = coverage(volumes, pairs)
    assert occupied[-1] <= 1 << 100
    numerator = sum((value - (occupied[t - 1] if t else 0)) * 9 ** (200 - t)
                    for t, value in enumerate(occupied))
    return numerator, occupied


def address_and_lane_checks():
    assert (0 < PARITY_BASE < TAG_BASE < GAMMA_BASE
            < PROFILE_BASE < PACKED_BASE < END <= LEADER_BASE == 1 << 31)
    assert len({displacement(t) for t in range(SLOTS)}) == SLOTS
    assert all(0 <= column < 1 << 31 for column in COLUMNS)
    gamma_values = [gamma(s) for s in range(SLOTS)]
    assert gamma_values[0] == 0
    for s,d in enumerate(gamma_values):
        assert gamma_values[d] == (s ^ d)
        if s:
            assert d not in (0,s)
    for bit in range(12):
        for s in range(SLOTS):
            assert gamma_values[s ^ (1<<bit)] == gamma_values[s] ^ gamma_values[1<<bit]
    scalar_addresses = {profile_address(i,t) for i in range(4) for t in range(SLOTS)}
    packed_addresses = {
        packed_address(i,p,g)
        for i in range(4) for p in range(8 if i else 1) for g in range(GROUPS)
    }
    assert scalar_addresses == set(range(PROFILE_BASE, PACKED_BASE))
    assert packed_addresses == set(range(PACKED_BASE,END))
    for b in range(126):
        for d in range(125):
            assert (((((d * REPLICATE) | HIGH) - b * REPLICATE) & HIGH) == HIGH) == (d >= b)
    rng = random.Random(1750)
    for _ in range(100):
        weights = [rng.randrange(125) for _ in range(8)]
        packed = sum(d << (8*i) for i,d in enumerate(weights))
        b = rng.randrange(126)
        assert (((packed | HIGH) - b * REPLICATE) & HIGH) == sum(
            128 << (8*i) for i,d in enumerate(weights) if d >= b
        )
    assert encode((1 << 1000) - 1).bit_length() <= 2000
    return len(scalar_addresses), len(packed_addresses)


def main():
    sector = json.load(sys.stdin)["sector"]
    hists = read_histograms()
    check_sector(sector, hists)
    small = {"sector":sector,"checked_tags":SECTOR_BOUNDARIES[sector+1]-SECTOR_BOUNDARIES[sector],"private_count":1 << PRIVATE_BITS}
    if sector == 3:
        small["scratch_lengths"] = address_and_lane_checks()
        numerator, occupied = component_numerator(hists)
        actual = 1 - Fraction(numerator ** BLOCKS, 10 ** 2000)
        manifest = json.loads((HERE / "submission.json").read_text())
        error_bound, log_bound = [
            Fraction(int(claim["score"]["numerator"]), int(claim["score"]["denominator"]))
            for claim in manifest["claims"]
        ]
        ops = operation_bound()
        assert actual <= error_bound and ops ** log_bound.denominator <= 2 ** log_bound.numerator
        small.update({
            "occupied_syndromes_lower_bounds": occupied,
            "component_success_numerator": str(numerator),
            "error_upper_bound": str(error_bound),"ops":ops,"log_ops_bound":str(log_bound),
        })
    print(json.dumps(small))


if __name__ == "__main__":
    main()
