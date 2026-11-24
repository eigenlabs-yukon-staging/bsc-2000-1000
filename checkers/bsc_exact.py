"""Organizer-owned exact-arithmetic checker for the BSC(1/10), n = 2000 challenge.

It computes one of two code-independent probabilities as an exact rational and
compares it with a requester-supplied rational bound B:

  tail(t)      = Pr[ Binomial(2000, 1/10) >= t ]
                 the chance that at least t of the 2000 bits flip.
  pairwise(w)  = Pr[ Binomial(w, 1/10) >= ceil(w/2) ]
                 an upper bound on the chance that ML decoding prefers one fixed
                 competing codeword at Hamming distance w, counting ties as errors
                 (so it holds for every tie-breaking rule).

A pass establishes only the arithmetic fact value <= B. It does not establish any
property of a submitted code (minimum distance, weight distribution, decoder).
"""
import json
import sys
from math import comb

N, P_NUM, P_DEN = 2000, 1, 10  # protected channel constants


def tail(n, t):
    # sum_{j >= t} C(n, j) * 1^j * 9^(n-j) / 10^n
    return sum(comb(n, j) * (P_DEN - P_NUM) ** (n - j) * P_NUM ** j for j in range(max(t, 0), n + 1)), P_DEN ** n


def decimal_up(num, den, digits=20):
    """Smallest decimal with `digits` significant digits that is >= num/den (0 <= num/den <= 1)."""
    if num == 0:
        return "0"
    e = 0
    while num * 10 ** (e + 1) < den * 10 ** digits:
        e += 1
    # scale so that the integer part has `digits` digits
    q, r = divmod(num * 10 ** e, den)
    if r:
        q += 1
    s = str(q).rjust(e + 1, "0")
    head, frac = s[:-e] if e else s, s[-e:] if e else ""
    return (head or "0") + ("." + frac if frac else "")


def parse_pos_int(s, name):
    if not (isinstance(s, str) and s.isascii() and s.isdigit() and (s == "0" or s[0] != "0")):
        raise ValueError(f"{name} must be a canonical decimal integer string")
    return int(s)


def main():
    d = json.load(sys.stdin)
    kind, k = d["kind"], d["parameter"]
    out = {"passed": False, "kind": kind, "parameter": k, "value_numerator": "",
           "value_denominator": "", "value_decimal_upper": "", "message": ""}
    try:
        bn = parse_pos_int(d["bound_numerator"], "bound_numerator")
        bd = parse_pos_int(d["bound_denominator"], "bound_denominator")
        if bd == 0:
            raise ValueError("bound_denominator must be positive")
        if kind == "tail":
            if not 0 <= k <= N + 1:
                raise ValueError("tail parameter t must satisfy 0 <= t <= 2001")
            num, den = tail(N, k)
        else:
            if not 1 <= k <= N:
                raise ValueError("pairwise parameter w must satisfy 1 <= w <= 2000")
            num, den = tail(k, (k + 1) // 2)
        out.update(value_numerator=str(num), value_denominator=f"10^{len(str(den)) - 1}",
                   value_decimal_upper=decimal_up(num, den), passed=num * bd <= bn * den,
                   message="exact value <= bound" if num * bd <= bn * den else "exact value exceeds bound")
    except (ValueError, KeyError) as exc:
        out["message"] = str(exc)
    print(json.dumps(out))


main()
