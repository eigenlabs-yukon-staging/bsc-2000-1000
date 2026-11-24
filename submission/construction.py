"""Explicit codes, information sets, table generation, and counted parameters."""
from hashlib import sha256
from itertools import combinations
from math import comb

COMPONENTS = ((127, 22, (4, 4, 4, 3), 15), (79, 14, (5, 4, 3), 1))
PAIRS = tuple(combinations(range(4), 2))


def polynomial_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def generator(p):
    residues = {x * x % p for x in range(1, p)}
    indicator = sum(1 << r for r in residues)
    polynomial = polynomial_gcd(indicator, (1 << p) | 1)
    extension = (polynomial.bit_count() & 1) << p
    rows = tuple((polynomial << i) | extension for i in range((p + 1) // 2))
    return polynomial, rows


def systematic(rows, coordinates):
    rows = list(rows)
    for i, coordinate in enumerate(coordinates):
        pivot = next(j for j in range(i, len(rows)) if (rows[j] >> coordinate) & 1)
        rows[i], rows[pivot] = rows[pivot], rows[i]
        for j in range(len(rows)):
            if j != i and (rows[j] >> coordinate) & 1:
                rows[j] ^= rows[i]
    return tuple(rows)


def standard_positions(p):
    k = (p + 1) // 2
    return tuple(range(k - 1, p))


def group_positions(p, seed):
    permutation = sorted(range(p + 1), key=lambda j: sha256(
        f"four-qr:{p}:{seed}:{j}".encode("ascii")
    ).digest())
    g = (p + 1) // 4
    return tuple(tuple(sorted(permutation[i:i + g])) for i in range(0, p + 1, g))


def paired_sets(groups):
    return tuple(groups[i] + groups[j] for i, j in PAIRS)


def linear_combination(rows, mask):
    word = 0
    while mask:
        bit = mask & -mask
        word ^= rows[bit.bit_length() - 1]
        mask ^= bit
    return word


def permitted(a, b, profile):
    small, large = min(a, b), max(a, b)
    return small < len(profile) and large <= profile[small]


def two_group_weights(g, profile):
    for a in range(g + 1):
        for b in range(g + 1):
            if permitted(a, b, profile):
                yield a, b


def table_entries(g, profile):
    return sum(comb(g, a) * comb(g, b) for a, b in two_group_weights(g, profile))


def delta_words(rows, profile):
    g = len(rows) // 2
    for a, b in two_group_weights(g, profile):
        for first in combinations(range(g), a):
            low_group = 0
            for i in first:
                low_group ^= rows[i]
            for second in combinations(range(g, 2 * g), b):
                word = low_group
                for i in second:
                    word ^= rows[i]
                yield word & ((1 << 64) - 1)
                yield word >> 64


def encode(message):
    sent = payload_offset = channel_offset = 0
    for p, _, _, repeat in COMPONENTS:
        n = p + 1
        k = n // 2
        _, rows = generator(p)
        basis = systematic(rows, standard_positions(p))
        for _ in range(repeat):
            mask = (message >> payload_offset) & ((1 << k) - 1)
            sent |= linear_combination(basis, mask) << channel_offset
            payload_offset += k
            channel_offset += n
    return sent


def sizes():
    words = 48
    operations = 0
    for p, _, profile, repeat in COMPONENTS:
        k = (p + 1) // 2
        scans = len(PAIRS)
        num = table_entries(k // 2, profile)
        words += scans * 2 * num
        operations += repeat * (2 + scans * (6 * k + 2 + 15 * num) + 5 + (k < 64) + 1)
    return words, operations
