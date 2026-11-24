import json, math, sys
from fractions import Fraction
json.load(sys.stdin)
exact_error=1-Fraction(9**1000,10**1000)
scale=10**100
success_floor=9**1000*scale//10**1000
B=1-Fraction(success_floor,scale)
assert exact_error<=B<Fraction(10**46-1,10**46)
assert math.gcd(B.numerator,B.denominator)==1
print(json.dumps({"length":2000,"message_bits":1000,"success_floor":str(success_floor),"B_numerator":str(B.numerator),"B_denominator":str(B.denominator),"exact_error_le_B":True,"fraction_is_reduced":True,"strictly_better_error_than_starter":True}))
