"""Fixed majority-plane schedule and exact arithmetic for its error."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import random
import sys

HERE = Path('/snapshot/files') if Path('/snapshot/files/opmodel.py').is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from opmodel import Machine

PLANES = 7


def encode(message):
    code = message
    for j in range(PLANES):
        plane = (message >> (64*j)) & (2**64-1)
        code |= plane << (64*(16+2*j))
        code |= plane << (64*(17+2*j))
    return code


def schedule():
    for j in range(16):
        yield 'load', 'T', j
        if j < PLANES:
            yield 'load', 'A', 16 + 2*j
            yield 'load', 'B', 17 + 2*j
            yield 'and_', 'C', 'T', 'A'
            yield 'or_', 'D', 'T', 'A'
            yield 'and_', 'D', 'D', 'B'
            yield 'or_', 'T', 'C', 'D'
        if j == 15:
            yield 'and_', 'T', 'T', 2**40-1
        yield 'store', None, 32 + j, 'T'


def run(received, program):
    ram = Machine([(received >> i) & 1 for i in range(2000)])
    registers = {}
    for op, dst, *args in program:
        values = (registers[x] if isinstance(x, str) else x for x in args)
        value = getattr(ram, op)(*values)
        if dst is not None:
            registers[dst] = value
    return sum(bit << i for i, bit in enumerate(ram.output_bits())), ram.ops


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / 'submission.json').read_text())
    b, ell = [Fraction(int(c['score']['numerator']), int(c['score']['denominator']))
              for c in manifest['claims']]
    p = Fraction(1, 10)
    correct = (1-p)**3 + 3*p*(1-p)**2
    success = correct**(64*PLANES)*(1-p)**(1000-64*PLANES)
    unit = 10**250
    assert b == 1 - Fraction(success.numerator*unit//success.denominator, unit)
    assert 1-success <= b < 1
    program = tuple(schedule())
    counts = Counter(op for op, *_ in program)
    w = sum(counts.values())
    assert w == 33+6*PLANES and w**ell.denominator <= 2**ell.numerator
    rng = random.Random(196)
    tests = [0, (1<<1000)-1] + [rng.getrandbits(1000) for _ in range(16)]
    for msg in tests:
        tx = encode(msg)
        assert run(tx, program) == (msg, w)
        y = tx ^ rng.getrandbits(2000)
        expected = (y & ((1<<1000)-1)) & ~((1<<(64*PLANES))-1)
        for j in range(64*PLANES):
            k, i = divmod(j, 64)
            copies = (j, 64*(16+2*k)+i, 64*(17+2*k)+i)
            expected |= (sum((y >> bit) & 1 for bit in copies) >= 2) << j
        assert run(y, program) == (expected, w)
    print(json.dumps({'bound': str(b), 'majority_success': str(correct),
                      'primitive_count': w, 'primitive_breakdown': counts, 'cases': 2*len(tests)}))


if __name__ == '__main__':
    main()
