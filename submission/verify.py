"""Exact distance, enumerator, error, and operation-certificate checks."""
from collections import Counter
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    LENGTH, ORDER, PRIME, SETS, WIDTH, bases, count_candidates, generator,
    golay_leaders, linear_combination, operation_bound, word_limit,
)


def packed_coordinates(word, positions):
    return sum(((word >> i) & 1) << j for j, i in enumerate(positions))


def verify_structure(rows, polynomial):
    assert polynomial.bit_length() == WIDTH
    assert all(packed_coordinates(row, SETS[0]) == 1 << i for i, row in enumerate(rows))
    assert all(row.bit_count() % 4 == 0 for row in rows)
    assert all((a & b).bit_count() % 2 == 0 for a in rows for b in rows)
    assert PRIME % 4 == 3 and all(PRIME % i != 0 for i in range(2, 11))
    translation = tuple((i + 1) % PRIME for i in range(PRIME)) + (PRIME,)
    inversion = (PRIME,) + tuple((-pow(i, -1, PRIME)) % PRIME for i in range(1, PRIME)) + (0,)
    for permutation in (translation, inversion):
        assert sorted(permutation) == list(range(LENGTH))
        for row in rows:
            moved = sum(((row >> i) & 1) << j for i, j in enumerate(permutation))
            mask = packed_coordinates(moved, SETS[0])
            assert linear_combination(rows, mask) == moved


def certify_distance(rows):
    """Search all systematic supports of sizes one through six."""
    def visit(start, remaining, word):
        if remaining == 1:
            for i in range(start, WIDTH):
                if (word ^ rows[i]).bit_count() <= 16:
                    return False
            return True
        for i in range(start, WIDTH - remaining + 1):
            if not visit(i + 1, remaining - 1, word ^ rows[i]):
                return False
        return True

    assert all(visit(0, weight, 0) for weight in range(1, 7))
    return sum(comb(WIDTH, weight) for weight in range(1, 7))


def krawtchouk(weight, at):
    return sum(
        (-1) ** i * comb(at, i) * comb(LENGTH - at, weight - i)
        for i in range(max(0, weight - (LENGTH - at)), min(weight, at) + 1)
    )


def solve_exact(matrix, variables):
    for i in range(variables):
        pivot = next(j for j in range(i, len(matrix)) if matrix[j][i] != 0)
        matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
        q = matrix[i][i]
        matrix[i] = [v / q for v in matrix[i]]
        for j, row in enumerate(matrix):
            if j != i and row[i]:
                q = row[i]
                matrix[j] = [x - q * y for x, y in zip(row, matrix[i])]
    assert all(all(v == 0 for v in row) for row in matrix[variables:])
    return [matrix[i][-1] for i in range(variables)]


def weight_distribution():
    # The all-ones codeword follows from self-duality and even weight.
    known = {0: 1, LENGTH: 1}
    unknowns = tuple(range(20, LENGTH - 16, 4))
    matrix = []
    for w in range(LENGTH + 1):
        coeff = [krawtchouk(w, d) - ((1 << WIDTH) if d == w else 0) for d in unknowns]
        rhs = (1 << WIDTH) * known.get(w, 0) - sum(
            c * krawtchouk(w, d) for d, c in known.items()
        )
        matrix.append([Fraction(v) for v in coeff + [rhs]])
    values = solve_exact(matrix, len(unknowns))
    assert all(v.denominator == 1 and v >= 0 for v in values)
    result = known | dict(zip(unknowns, map(int, values)))
    assert sum(result.values()) == 1 << WIDTH
    return result


def binomial_ball(n, t):
    return sum(comb(n, i) for i in range(min(n, t) + 1))


def full_intersection(distance, t):
    return sum(
        comb(distance, i) * binomial_ball(LENGTH - distance, t - max(i, distance - i))
        for i in range(distance + 1)
        if t >= max(i, distance - i)
    )


def restricted_volume(t):
    return sum(
        comb(WIDTH, i) * comb(WIDTH, j)
        for i in range(min(WIDTH, t) + 1)
        for j in range(min(WIDTH, t - i) + 1)
        if min(i, j) <= ORDER
    )


def component_success(distribution):
    previous = 0
    numerator = 0
    for t in range(23):
        volume = restricted_volume(t)
        collisions = sum(
            count * full_intersection(d, t)
            for d, count in distribution.items() if 0 < d <= 2 * t
        )
        denom = volume + collisions
        bonferroni = volume - collisions // 2
        cauchy = (volume * volume + denom - 1) // denom
        d = 1 + collisions // volume
        scale = d * (d + 1)
        secant = (2 * d * volume - collisions + scale - 1) // scale
        coverage = max(previous, bonferroni, cauchy, secant)
        assert 0 <= coverage <= 1 << WIDTH
        numerator += (coverage - previous) * 9 ** (LENGTH - t)
        previous = coverage
    return numerator


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    scores = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"])) for c in manifest["claims"]]
    polynomial, _ = generator()
    qr_rows, other_rows = bases()
    verify_structure(qr_rows, polynomial)
    assert all(packed_coordinates(row, SETS[1]) == 1 << i for i, row in enumerate(other_rows))
    searched = certify_distance(qr_rows)
    dist = weight_distribution()
    success = component_success(dist)
    tail_hist = Counter(v.bit_count() for v in golay_leaders().values())
    assert tail_hist == {0: 1, 1: 24, 2: 276, 3: 2024, 4: 1771}
    tail_success = sum(count * 9 ** (24 - w) for w, count in tail_hist.items())
    bound, exponent = scores
    assert 1 - Fraction(success ** 19 * tail_success, 10 ** 2000) <= bound
    worst = operation_bound()
    assert word_limit() <= 1 << 32
    assert worst ** exponent.denominator <= 2 ** exponent.numerator
    print(json.dumps({
        "qr_polynomial_hex": hex(polynomial), "systematic_candidates_checked": searched,
        "nonzero_weight_counts": dist, "code_distance_lower_bound": 20,
        "qr_success_numerator_lower_bound": success, "golay_histogram": tail_hist,
        "golay_success_numerator": tail_success, "block_error_upper_bound": str(bound),
        "candidates_per_information_set": count_candidates(), "ops_upper_bound": worst,
        "memory_word_limit": word_limit(), "log_ops_upper_bound": str(exponent),
    }))


if __name__ == "__main__":
    main()
