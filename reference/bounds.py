"""Organizer reference: exact anchors for a (2000, 2^1000) code on BSC(1/10).

Prints the sphere-packing converse (a lower bound on the average block error
probability of EVERY code with these parameters, under any decoder) and the
random-coding union (RCU) bound (an upper bound on the error averaged over the
ensemble of i.i.d. uniform codebooks with ML decoding, ties counted as errors).
All arithmetic is exact (integers and fractions); only the printed decimals
are rounded. Python 3 standard library only. Run: python3 -I reference/bounds.py
"""
from fractions import Fraction
from math import comb

N, K = 2000, 1000
# P(Z has weight w) = C(N,w) * 9^(N-w) / 10^N  for crossover 1/10
DEN = 10 ** N
BIN = [comb(N, w) for w in range(N + 1)]
PW = [BIN[w] * 9 ** (N - w) for w in range(N + 1)]


def converse():
    # Correct-decoding probability of a message whose decoding region has
    # volume v is at most f(v): the probability of the v most likely error
    # patterns (lowest weight first, fractional on the boundary shell).
    # f is concave and increasing; regions are disjoint, so the average
    # correct probability is at most f(2^N / M) = f(2^(N-K)).
    vol = 2 ** (N - K)
    used, r = 0, 0
    while used + BIN[r] <= vol:
        used += BIN[r]
        r += 1
    frac = Fraction(vol - used, BIN[r])
    correct = Fraction(sum(PW[:r]), DEN) + frac * Fraction(PW[r], DEN)
    return 1 - correct, r, frac


def rcu():
    m = 2 ** K
    acc, cum = Fraction(0), 0
    for t in range(N + 1):
        cum += BIN[t]
        acc += Fraction(PW[t], DEN) * min(Fraction(1), Fraction((m - 1) * cum, 2 ** N))
    return acc


SCORE_MINIMUM = Fraction(525481183, 10 ** 10)

if __name__ == "__main__":
    conv, r, frac = converse()
    print(f"converse: every code has error >= {float(conv):.12f} (shell radius {r}, boundary fraction {float(frac):.6f})")
    print(f"protected score_minimum {SCORE_MINIMUM} = {float(SCORE_MINIMUM):.12f}; above converse: {SCORE_MINIMUM > conv}")
    print(f"random-coding (RCU, ensemble) error <= {float(rcu()):.12f}")
