"""Explicit QR(103) basis, systematic information sets, and decoder sizes."""
from itertools import combinations
from math import comb

PRIME = 103
LENGTH = PRIME + 1
WIDTH = LENGTH // 2
ORDER = 6
SETS = (tuple(range(51, 103)), tuple(range(51)) + (103,))
GOLAY_POLY = sum(1 << i for i in (0, 2, 4, 5, 6, 10, 11))


def polynomial_gcd(a, b):
    while b:
        while a.bit_length() >= b.bit_length():
            a ^= b << (a.bit_length() - b.bit_length())
        a, b = b, a
    return a


def systematic(rows, coordinates):
    rows = list(rows)
    for i, coordinate in enumerate(coordinates):
        pivot = next(j for j in range(i, len(rows)) if (rows[j] >> coordinate) & 1)
        rows[i], rows[pivot] = rows[pivot], rows[i]
        for j in range(len(rows)):
            if j != i and (rows[j] >> coordinate) & 1:
                rows[j] ^= rows[i]
    return tuple(rows)


def generator():
    residues = {x * x % PRIME for x in range(1, PRIME)}
    indicator = sum(1 << r for r in residues)
    poly = polynomial_gcd(indicator, (1 << PRIME) | 1)
    rows = tuple((poly << i) | ((poly.bit_count() & 1) << PRIME) for i in range(WIDTH))
    return poly, rows


def bases():
    _, rows = generator()
    return tuple(systematic(rows, positions) for positions in SETS)


def linear_combination(rows, mask):
    word = 0
    while mask:
        bit = mask & -mask
        word ^= rows[bit.bit_length() - 1]
        mask ^= bit
    return word


def qr_delta_words(rows):
    for weight in range(ORDER + 1):
        for indices in combinations(range(WIDTH), weight):
            word = 0
            for i in indices:
                word ^= rows[i]
            yield word & ((1 << 64) - 1)
            yield word >> 64


def golay_encode(data):
    encoded = linear_combination(tuple(GOLAY_POLY << i for i in range(12)), data)
    return encoded | ((encoded.bit_count() & 1) << 23)


def golay_syndrome(word):
    par = word.bit_count() & 1
    body = word & ((1 << 23) - 1)
    for i in range(22, 10, -1):
        if (body >> i) & 1:
            body ^= GOLAY_POLY << (i - 11)
    return body | (par << 11)


def golay_leaders():
    leaders = {}
    for weight in range(5):
        for positions in combinations(range(24), weight):
            word = sum(1 << i for i in positions)
            leaders.setdefault(golay_syndrome(word), word)
    return leaders


def golay_table_words():
    leaders = golay_leaders()
    table = [0] * (1 << 24)
    for data in range(1 << 12):
        encoded = golay_encode(data)
        for noise in leaders.values():
            table[encoded ^ noise] = data
    return table


def encode(message):
    qr_rows, _ = bases()
    out = 0
    for i in range(19):
        data = (message >> (WIDTH * i)) & ((1 << WIDTH) - 1)
        out |= linear_combination(qr_rows, data) << (LENGTH * i)
    return out | (golay_encode(message >> (19 * WIDTH)) << (19 * LENGTH))


def count_candidates():
    return sum(comb(WIDTH, i) for i in range(ORDER + 1))


def word_limit():
    return 48 + 4 * count_candidates() + (1 << 24)


def operation_bound():
    return 19 * (649 + 30 * count_candidates()) + 27
