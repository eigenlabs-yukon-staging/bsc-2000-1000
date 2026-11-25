"""Exact arithmetic for the two declared scores."""
import json
import math
import sys


def rational_input(values):
    numerator, denominator = map(int, values)
    assert denominator > 0 and math.gcd(numerator, denominator) == 1
    return numerator, denominator


data = json.load(sys.stdin)
b, d = rational_input(data["B"])
p, q = rational_input(data["L"])
error_passed = (d - b) * 2**808 * 10**192 <= d * 9**192
operations_passed = 19**q <= 2**p
assert error_passed and operations_passed
print(json.dumps({"error_passed": error_passed,
                  "operations_passed": operations_passed,
                  "B": data["B"], "L": data["L"]}))
