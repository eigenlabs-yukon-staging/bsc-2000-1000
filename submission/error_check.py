"""Count the syndromes represented at each certified noise weight."""
import json
from fractions import Fraction
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import syndrome_columns


def histogram(columns, width, low, radius):
    """Bitset of reachable syndromes, as rows indexed by their high bits."""
    slots = 1 << (width - low)
    mask = (1 << low) - 1
    row_bits = 1 << low
    swap_masks = []
    for j in range(low):
        step = 1 << j
        lower_lanes = ((1 << row_bits) - 1) // ((1 << step) + 1)
        swap_masks.append((step, lower_lanes, lower_lanes << step))
    swaps = [
        tuple(swap_masks[j] for j in range(low) if (value >> j) & 1)
        for value in range(1 << low)
    ]
    visited = [0] * slots
    visited[0] = 1
    previous = 1
    observed = [1]
    edge_order = sorted(columns, key=lambda col: (col & mask).bit_count())
    full = (1 << row_bits) - 1
    for _ in range(radius):
        reached = visited.copy()
        if 4 * previous >= 3 * (1 << width):
            for i, union in enumerate(visited):
                if union != full:
                    for col in edge_order:
                        row = visited[i ^ (col >> low)]
                        for shift, lo, hi in swaps[col & mask]:
                            row = ((row & lo) << shift) | ((row & hi) >> shift)
                        union |= row
                        if union == full:
                            break
                reached[i] = union
        else:
            expand_columns(columns, low, mask, swaps, visited, reached)
        count = sum(row.bit_count() for row in reached)
        observed.append(count - previous)
        if count == 1 << width:
            break
        previous = count
        visited = reached
    return observed


def expand_columns(columns, low, mask, swaps, visited, reached):
    for col in columns:
        perm = swaps[col & mask]
        high = col >> low
        if perm:
            for i, row in enumerate(visited):
                if row:
                    for shift, lo, hi in perm:
                        row = ((row & lo) << shift) | ((row & hi) >> shift)
                    reached[i ^ high] |= row
        else:
            for i, row in enumerate(visited):
                reached[i ^ high] |= row

def weighted_count(hist, k):
    return sum(count * 9 ** (2 * k - weight) for weight, count in enumerate(hist))


def main():
    json.load(sys.stdin)
    claim = json.loads((HERE / "submission.json").read_text())["claims"][0]["score"]
    bound = Fraction(int(claim["numerator"]), int(claim["denominator"]))
    wide = histogram(syndrome_columns(30), 30, 13, 9)
    tail = histogram(syndrome_columns(10), 10, 5, 10)
    assert sum(tail) == 1 << 10
    a = weighted_count(wide, 30)
    b = weighted_count(tail, 10)
    exact_upper_bound = 1 - Fraction(a ** 33 * b, 10 ** 2000)
    assert exact_upper_bound <= bound
    print(json.dumps({"big_histogram": wide, "tail_histogram": tail,
                      "big_weighted_count": a, "tail_weighted_count": b,
                      "error_upper_bound": str(bound), "bound_verified": True}))


if __name__ == "__main__":
    main()
