import json
import sys

d = json.load(sys.stdin)["exponent"]
lhs = 9 ** 1000 * 10 ** d
rhs = 10 ** 1000
print(json.dumps({"exponent": d, "nine_pow_1000_times_ten_pow_exponent_ge_ten_pow_1000": lhs >= rhs,
                  "lhs_decimal_digits": len(str(lhs)), "rhs_decimal_digits": len(str(rhs))}))
