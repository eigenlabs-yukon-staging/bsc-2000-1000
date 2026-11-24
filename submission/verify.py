"""Exact calculation and replay of the nineteen-primitive decoder."""
from fractions import Fraction
import json
from pathlib import Path
import random
import sys


HERE = Path('/snapshot/files') if Path('/snapshot/files/opmodel.py').is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from opmodel import Machine


def encode(message):
    return message | (message << 1000)


def decode(received):
    ram = Machine([(received >> i) & 1 for i in range(2000)])
    for i in range(3):
        ram.store(32 + i, ram.load(i))
    for i in range(3, 16):
        ram.store(32 + i, 0)
    return sum(bit << i for i, bit in enumerate(ram.output_bits())), ram.ops


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / 'submission.json').read_text())
    bound, exponent = [Fraction(int(c['score']['numerator']), int(c['score']['denominator']))
                       for c in manifest['claims']]
    success = Fraction(9, 10)**192 * Fraction(1, 2)**808
    unit = 10**255
    assert bound == 1 - Fraction(success.numerator*unit//success.denominator, unit)
    assert 1-success <= bound < 1
    rng = random.Random(613)
    cases = [0, 2**2000 - 1]
    cases.extend(encode(rng.getrandbits(1000)) ^ rng.getrandbits(2000) for _ in range(16))
    for received in cases:
        output, ops = decode(received)
        assert output == received % 2**192
        assert ops == 19 and ops**exponent.denominator <= 2**exponent.numerator
    print(json.dumps({'success': str(success), 'rounded_error': str(bound),
                      'operations': ops, 'replayed_words': len(cases)}))


if __name__ == '__main__':
    main()
