"""Exact syndrome histograms with short-mask and bitset expansions."""
from fractions import Fraction
from functools import reduce
from itertools import combinations
import json
from operator import xor
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import syndrome_columns

BASIS_ROWS = (36, 80, 128, 160, 256, 576, 803, 1408, 1417, 2083, 2305, 2834, 4320, 8226, 16389, 32962, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912)
EDGE_LIMITS = (28, 32, 36, 40, 60)


def check_basis():
    pivots = {}
    for initial in BASIS_ROWS:
        assert 0 < initial < 1 << 30
        row = initial
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    assert len(BASIS_ROWS) == len(pivots) == 30


def expand(rows, new, selected, permutations, low, full, count):
    mask = (1 << low) - 1
    edges = [(c >> low, permutations[c & mask]) for c in selected]
    if 4 * count >= 3 * (1 << 30):
        for j, union in enumerate(rows):
            if union == full:
                continue
            for high, swaps in edges:
                row = rows[j ^ high]
                for shift, left, right in swaps:
                    row = ((row & left) << shift) | ((row & right) >> shift)
                union |= row
                if union == full:
                    break
            new[j] = union
    else:
        for high, swaps in edges:
            if not swaps:
                for j, row in enumerate(rows):
                    new[j ^ high] |= row
                continue
            for j, row in enumerate(rows):
                if row:
                    for shift, left, right in swaps:
                        row = ((row & left) << shift) | ((row & right) >> shift)
                    new[j ^ high] |= row


def mark_combinations(rows, columns, weight, low, singleton):
    mask = (1 << low) - 1
    for group in combinations(columns, weight):
        syndrome = reduce(xor, group, 0)
        rows[syndrome >> low] |= singleton[syndrome & mask]


def wide_histogram():
    check_basis()
    columns = tuple(sum(((value & row).bit_count() & 1) << j
                        for j, row in enumerate(BASIS_ROWS))
                    for value in syndrome_columns(30))
    low = 13
    mask = (1 << low) - 1
    full = (1 << (1 << low)) - 1
    swap_masks = []
    for j in range(low):
        shift = 1 << j
        left = full // ((1 << shift) + 1)
        swap_masks.append((shift, left, left << shift))
    permutations = [tuple(swap_masks[j] for j in range(low) if (v >> j) & 1)
                    for v in range(1 << low)]
    order = sorted(columns, key=lambda c: (c & mask).bit_count())
    singleton = tuple(1 << j for j in range(1 << low))
    rows = [0] * (1 << (30-low))
    previous = 0
    histogram = []
    for weight in range(6):
        mark_combinations(rows, columns, weight, low, singleton)
        count = sum(row.bit_count() for row in rows)
        histogram.append(count - previous)
        previous = count
    for weight, limit in enumerate(EDGE_LIMITS, 6):
        new = rows.copy()
        expand(rows, new, order[:limit], permutations, low, full, previous)
        mark_combinations(new, order[limit:], weight, low, singleton)
        count = sum(row.bit_count() for row in new)
        histogram.append(count - previous)
        previous = count
        rows = new
        if count == 1 << 30:
            break
    assert previous == 1 << 30
    return histogram


def tail_histogram():
    visited = set()
    columns = syndrome_columns(10)
    result = []
    for weight in range(11):
        count = len(visited)
        visited.update(reduce(xor, group, 0) for group in combinations(columns, weight))
        result.append(len(visited) - count)
        if len(visited) == 1 << 10:
            break
    assert len(visited) == 1 << 10
    return result


def weighted_count(histogram, width):
    return sum(count * 9 ** (2 * width - weight)
               for weight, count in enumerate(histogram))


def main():
    json.load(sys.stdin)
    claim = json.loads((HERE / "submission.json").read_text())["claims"][0]["score"]
    bound = Fraction(int(claim["numerator"]), int(claim["denominator"]))
    wide, tail = wide_histogram(), tail_histogram()
    assert wide == [1, 60, 1770, 34220, 487635, 5460504, 49216586, 321043162, 647670432, 49827096, 358]
    assert tail == [1, 20, 190, 740, 71, 2]
    a, b = weighted_count(wide, 30), weighted_count(tail, 10)
    error = 1 - Fraction(a ** 33 * b, 10 ** 2000)
    assert error <= bound
    print(json.dumps({"wide_histogram": wide, "tail_histogram": tail, "basis_rank": 30,
                      "wide_weighted_count": a, "tail_weighted_count": b,
                      "error_upper_bound": str(bound), "bound_verified": True}))


if __name__ == "__main__":
    main()
