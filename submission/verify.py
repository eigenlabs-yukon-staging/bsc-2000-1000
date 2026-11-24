"""Pinned exact searches and algebra for the two explicit components."""
from collections import deque
from fractions import Fraction
import json
from math import comb, isqrt
from pathlib import Path
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    COMPONENTS, generator, group_positions, linear_combination, paired_sets,
    permitted, sizes, standard_positions, systematic, table_entries,
)


def packed(word, positions):
    return sum(((word >> i) & 1) << j for j, i in enumerate(positions))


def automorphisms(p):
    assert p % 4 == 3 and all(p % i != 0 for i in range(2, isqrt(p) + 1))
    translation = tuple((i + 1) % p for i in range(p)) + (p,)
    inversion = (p,) + tuple((-pow(i, -1, p)) % p for i in range(1, p)) + (0,)
    return translation, inversion


def verify_structure(p, rows, poly, seed):
    n = p + 1
    k = n // 2
    assert poly.bit_length() == k
    assert all(row.bit_count() % 4 == 0 for row in rows)
    assert all((a & b).bit_count() % 2 == 0 for a in rows for b in rows)
    standard = standard_positions(p)
    basis = systematic(rows, standard)
    assert all(packed(row, standard) == 1 << i for i, row in enumerate(basis))
    groups = group_positions(p, seed)
    assert len(groups) == 4 and all(len(group) == n // 4 for group in groups)
    assert sorted(j for group in groups for j in group) == list(range(n))
    for positions in paired_sets(groups):
        assert all(packed(row, positions) == 1 << i
                   for i, row in enumerate(systematic(rows, positions)))
    for perm in automorphisms(p):
        assert sorted(perm) == list(range(n))
        for row in basis:
            moved = sum(((row >> i) & 1) << j for i, j in enumerate(perm))
            assert linear_combination(basis, packed(moved, standard)) == moved
    return basis


def scan_basis(basis, degree, lower, target=None):
    k = len(basis)
    hits = []
    def visit(start, remaining, word):
        if remaining == 1:
            for i in range(start, k):
                final = word ^ basis[i]
                weight = final.bit_count()
                if weight < lower:
                    raise ValueError("unexpected low-weight codeword")
                if weight == target:
                    hits.append(final)
            return
        for i in range(start, k - remaining + 1):
            visit(i + 1, remaining - 1, word ^ basis[i])
    for w in range(1, degree + 1):
        visit(0, w, 0)
    return hits, sum(comb(k, w) for w in range(1, degree + 1))


def orbit_count(p, hits):
    _, inversion = automorphisms(p)
    targets = tuple(1 << j for j in inversion)
    top = 1 << p
    seen = set()
    orbits = []
    for root in hits:
        if root in seen:
            continue
        before = len(seen)
        seen.add(root)
        queue = deque([root])
        while queue:
            word = queue.popleft()
            rotated = ((word << 1) & (top - 1)) | ((word >> (p - 1)) & 1) | (word & top)
            reciprocal = 0
            z = word
            while z:
                bit = z & -z
                reciprocal |= targets[bit.bit_length() - 1]
                z ^= bit
            for moved in (rotated, reciprocal):
                if moved not in seen:
                    seen.add(moved)
                    queue.append(moved)
        orbits.append(len(seen) - before)
    return len(seen), orbits


def krawtchouk(n, weight, at):
    return sum((-1) ** i * comb(at, i) * comb(n - at, weight - i)
               for i in range(max(0, weight - n + at), min(weight, at) + 1))


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


def weight_distribution(n, lowest, extras):
    k = n // 2
    known = {0: 1, n: 1}
    for w, count in extras.items():
        known[w] = known[n - w] = count
    unknowns = [w for w in range(lowest, n - lowest + 1, 4) if w not in known]
    matrix = []
    for w in range(n + 1):
        coeff = [krawtchouk(n, w, d) - ((1 << k) if d == w else 0) for d in unknowns]
        rhs = (1 << k) * known.get(w, 0) - sum(
            c * krawtchouk(n, w, d) for d, c in known.items()
        )
        matrix.append([Fraction(v) for v in coeff + [rhs]])
    values = solve_exact(matrix, len(unknowns))
    assert all(v.denominator == 1 and v >= 0 for v in values)
    result = known | dict(zip(unknowns, map(int, values)))
    assert sum(result.values()) == 1 << k
    return result


def ball(n, t):
    return sum(comb(n, i) for i in range(min(n, t) + 1))


def intersection(n, d, t):
    return sum(comb(d, i) * ball(n - d, t - max(i, d - i))
               for i in range(max(0, d - t), min(d, t) + 1))


def component_success(n, distribution, profile, tmax):
    g = n // 4
    choose = [comb(g, w) for w in range(min(g, tmax) + 1)]
    shells = [0] * (tmax + 1)
    for a in range(len(choose)):
        for b in range(min(g, tmax - a) + 1):
            for c in range(min(g, tmax - a - b) + 1):
                ab = choose[a] * choose[b] * choose[c]
                for d in range(min(g, tmax - a - b - c) + 1):
                    first, second = sorted((a, b, c, d))[:2]
                    if permitted(first, second, profile):
                        shells[a + b + c + d] += ab * choose[d]
    previous = volume = output = 0
    for t, shell in enumerate(shells):
        volume += shell
        m = sum(count * intersection(n, d, t)
                for d, count in distribution.items() if 0 < d <= 2 * t)
        secant = 1 + m // volume
        denom = secant * (secant + 1)
        lower = max(previous, volume - m // 2, -(-(volume * volume) // (volume + m)),
                    -(-(2 * secant * volume - m) // denom))
        assert previous <= lower <= 1 << (n // 2)
        output += (lower - previous) * 9 ** (n - t)
        previous = lower
    return output


def verify_component(p, seed, profile):
    poly, rows = generator(p)
    basis = verify_structure(p, rows, poly, seed)
    if p == 127:
        distance = 20
        hits, searched = scan_basis(basis, 7, distance, 20)
        coefficient, orbits = orbit_count(p, hits)
        known = {20: coefficient}
    else:
        distance = 16
        hits, searched = scan_basis(basis, 4, distance)
        known, orbits = {}, []
    n = p + 1
    distribution = weight_distribution(n, distance, known)
    success = component_success(n, distribution, profile, 22 if p == 127 else 18)
    return success, {"polynomial_hex": hex(poly), "searched": searched, "hits": len(hits),
                     "orbit_sizes": orbits, "known_low_counts": known,
                     "weight_distribution": distribution, "success_numerator": success,
                     "candidates_per_scan": table_entries(n // 4, profile)}


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    B, L = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
            for c in manifest["claims"]]
    numerator = 1
    results = {}
    for p, seed, profile, repeated in COMPONENTS:
        success, data = verify_component(p, seed, profile)
        numerator *= success ** repeated
        results[p] = data
    assert 1 - Fraction(numerator, 10 ** 2000) <= B
    words, operations = sizes()
    assert words <= 1 << 32
    assert operations ** L.denominator <= 2 ** L.numerator
    print(json.dumps({"components": results, "error_bound": str(B), "memory_words": words,
                      "ops": operations, "log_ops_bound": str(L)}))


if __name__ == "__main__":
    main()
