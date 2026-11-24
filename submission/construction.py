"""The explicit four-group field code and its fixed tables."""
from functools import lru_cache
from itertools import combinations
from math import comb

M = 500
MASK = (1 << M) - 1
MOD_LOW = 1 | (1 << 8) | (1 << 13) | (1 << 21)
MODULUS = (1 << M) | MOD_LOW
LAMBDA = 3
THIRD = MASK // 3
CUBE_EXP = (2 * THIRD + 1) // 3
RADIUS = 108
BYTES = 63
PAIRS = tuple(combinations(range(4), 2))
CHART_COSTS = (30000, 165000, 165000, 165000, 165000, 725000)
ARRAY_LIST = (
    "Y0", "Y1", "Y2", "Y3", "KEY0", "KEY1", "V0", "V1", "V2", "V3",
    "BEST0", "BEST1", "X2", "X4", "T", "G", "H", "J", "CHECK", "C1", "C2",
    "INV", "INV2", "TEMP", "K", "TARGET", "ASROOT", "A0", "A1", "BASEA", "BASEB",
    "ROOT", "INVA",
)
BASE = {name: 48 + 8 * i for i, name in enumerate(ARRAY_LIST)}
NEXT = 48 + 8 * len(ARRAY_LIST)
LAST = NEXT + RADIUS + 1
MULT = LAST + RADIUS + 1
REDSHIFT1 = MULT + 16 * 8
REDSHIFT4 = REDSHIFT1 + 2
MAP_BASE = REDSHIFT4 + 16


def remainder(a, b=MODULUS):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def carryless(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def product(a, b):
    return remainder(carryless(a, b))


def square(a):
    return product(a, a)


def cube(a):
    return product(a, square(a))


def frobenius(a, s):
    for _ in range(s):
        a = square(a)
    return a


def power(a, e):
    result = 1
    for bit in bin(e)[2:]:
        result = square(result)
        if bit == "1":
            result = product(result, a)
    return result


def geom(a, length, step):
    result = a
    k = 1
    for bit in bin(length)[3:]:
        result = product(result, frobenius(result, step * k))
        if bit == "1":
            result = product(result, frobenius(a, 2 * step * k))
        k = 2 * k + int(bit)
    return result


def inverse(a):
    assert a
    return square(geom(a, M - 1, 1))


def cube_root(a):
    a7 = product(product(square(a), frobenius(a, 2)), a)
    t = frobenius(a7, 3)
    v = product(a, geom(t, 83, 6))
    return v if cube(v) == a else None


@lru_cache(None)
def constant(tag):
    if tag == "lambda":
        return LAMBDA
    if tag == "lambda_inv":
        return inverse(LAMBDA)
    if tag == "zeta":
        return power(LAMBDA, THIRD)
    if tag == "zeta2":
        return square(constant("zeta"))
    raise ValueError("unknown field constant")


@lru_cache(None)
def AS_pivots():
    pivots = {}
    for i in range(M):
        preimage = 1 << i
        image = square(preimage) ^ preimage
        while image:
            pivot = image.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (image, preimage)
                break
            r, x = pivots[pivot]
            image ^= r
            preimage ^= x
    return tuple((i, *pivots[i]) for i in sorted(pivots, reverse=True))


def AS_map(y):
    x = 0
    for pivot, image, preimage in AS_pivots():
        if (y >> pivot) & 1:
            y ^= image
            x ^= preimage
    return x


def cube_roots(x):
    base = cube_root(x)
    if base is None:
        return ()
    return (base, product(constant("zeta"), base), product(constant("zeta2"), base))


def codeword(a, b):
    return a, b, product(a, b), cube(a) ^ product(LAMBDA, cube(b))


def encode(message):
    coordinates = codeword(message & MASK, message >> M)
    return sum(value << (M * i) for i, value in enumerate(coordinates))


def chart(pair, x, y):
    """Host semantic counterpart of the bounded field chart macros."""
    if pair == (0, 1):
        return ((x, y),) if x and y else ()
    if not x:
        return ()
    if pair == (0, 2):
        return ((x, product(y, inverse(x))),) if y else ()
    if pair == (1, 2):
        return ((product(y, inverse(x)), x),) if y else ()
    if pair == (0, 3):
        target = product(y ^ cube(x), constant("lambda_inv"))
        return tuple((x, b) for b in cube_roots(target)) if target else ()
    if pair == (1, 3):
        target = y ^ product(LAMBDA, cube(x))
        return tuple((a, x) for a in cube_roots(target)) if target else ()
    assert pair == (2, 3)
    c, d = x, y
    k = product(LAMBDA, cube(c))
    if not d:
        possible = (frobenius(k, M - 1),)
    else:
        z = product(k, square(inverse(d)))
        w = AS_map(z)
        if square(w) ^ w != z:
            return ()
        first = product(d, w)
        possible = (first, first ^ d)
    result = []
    for value in possible:
        a = cube_root(value)
        if a is None:
            continue
        b = product(c, inverse(a))
        result.append((a, b))
        result.append((product(constant("zeta"), a), product(constant("zeta2"), b)))
        result.append((product(constant("zeta2"), a), product(constant("zeta"), b)))
    return tuple(result)


def map_tags():
    shifts = {1, 2, 3, M - 1}
    for length, step in ((83, 6), (M - 1, 1)):
        k = 1
        for bit in bin(length)[3:]:
            shifts.add(step * k)
            if bit == "1":
                shifts.add(2 * step * k)
            k = 2 * k + int(bit)
    return tuple(sorted(shifts)) + ("lambda", "lambda_inv", "zeta", "zeta2", "AS")


def map_columns(tag):
    if tag == "AS":
        for i in range(M):
            yield AS_map(1 << i)
    else:
        power_base = frobenius(2, tag) if isinstance(tag, int) else 2
        value = 1 if isinstance(tag, int) else constant(tag)
        for _ in range(M):
            yield value
            value = product(value, power_base)
    for _ in range(8 * BYTES - M):
        yield 0


def map_words(tag):
    cols = iter(map_columns(tag))
    for _ in range(BYTES):
        vals = [0]
        for _ in range(8):
            v = next(cols)
            vals.extend(x ^ v for x in vals[:])
        for value in vals:
            for h in range(8):
                yield (value >> (64 * h)) & ((1 << 64) - 1)


def reduction_table(s):
    return tuple(carryless(c, MOD_LOW) for c in range(1 << s))


MEMORY_END = MAP_BASE + len(map_tags()) * BYTES * 256 * 8


def W_bound():
    n = sum(comb(2 * M, i) for i in range(RADIUS + 1))
    return n, (sum(CHART_COSTS) + 6 * 48) * n + 600
