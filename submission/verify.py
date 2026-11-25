"""Verify the two exact score inequalities for the affine-tail code."""
from collections import Counter
from itertools import combinations
import json
import math
import sys


def rational_input(values):
    numerator, denominator = map(int, values)
    assert denominator > 0 and math.gcd(numerator, denominator) == 1
    return numerator, denominator


columns = [*range(1, 16), *range(16, 256, 16), 17, 34]
weights = {0: 0}
for c in columns:
    weights.setdefault(c, 1)
for a, b in combinations(columns, 2):
    weights.setdefault(a ^ b, 2)
assert len(weights) == 256
hist = Counter(weights.values())
assert hist == {0: 1, 1: 32, 2: 223}

excluded = Counter()
for v in range(48):
    s = 248
    for i, c in enumerate(columns):
        if (v >> i) & 1:
            s ^= c
    excluded[weights[s]] += 1
assert excluded == {2: 48}

whole = sum(count * 9**(32 - w) for w, count in hist.items())
lost = sum(count * 9**(32 - w) for w, count in excluded.items())
success_num = ((whole << 24) - lost) * 9**976
success_den = (1 << 24) * 10**1008
data = json.load(sys.stdin)
b, d = rational_input(data["B"])
p, q = rational_input(data["L"])
error_passed = (d - b) * success_den <= d * success_num
operations_passed = 37**q <= 2**p
assert error_passed and operations_passed
print(json.dumps({"error_passed": error_passed,
                  "operations_passed": operations_passed,
                  "weights": hist, "excluded_weights": excluded,
                  "B": data["B"], "L": data["L"]}))
