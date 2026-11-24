"""Verify excluded-address probability and evaluate the fixed RAM schedule."""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import json
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    MAIN_COLUMNS, TAIL_COLUMNS, TAIL_BASE, encode, leaders, parity_table, schedule,
    table_word,
)
from opmodel import Machine


def success_counts(main_parity, main_leaders, tail_leaders):
    histogram = Counter(e.bit_count() for e in main_leaders)
    tail_histogram = Counter(e.bit_count() for e in tail_leaders)
    numerator = sum(c * 9 ** (32 - w) for w, c in histogram.items())
    # For a successful main block, a given leader e and a given parity
    # half of the received word determine the target transmitted parity.
    multiplicities = Counter(main_parity)
    stripe = sum(
        multiplicities[65535 ^ (e >> 16)] * 9 ** (32 - e.bit_count())
        for e in main_leaders
    )
    prefix = sum(9 ** (32 - main_leaders[main_parity[r]].bit_count()) for r in range(48))
    safe_main = Fraction((1 << 16) * numerator - stripe - prefix, (1 << 16) * 10 ** 32)
    tail = Fraction(sum(c * 9 ** (16 - w) for w, c in tail_histogram.items()), 10 ** 16)
    return safe_main, tail, histogram, tail_histogram, stripe, prefix


class SampleMachine(Machine):
    def load(self, address):
        if address >= 48 and address not in self.mem:
            raise ValueError("sample accessed an entry that was not preloaded")
        return super().load(address)


def replay(machine, program):
    regs = {}
    for op, dst, *operands in program:
        args = [regs[v] if isinstance(v, str) else v for v in operands]
        if op == "move":
            regs[dst] = args[0]
        elif op == "store":
            machine.store(*args)
        else:
            regs[dst] = getattr(machine, op)(*args)
    return machine.output_bits()


def bits(word, width):
    return [(word >> i) & 1 for i in range(width)]


def sample(received, program, table_args, cost):
    keys = [(received >> (32 * j)) & 0xFFFFFFFF for j in range(62)]
    tail_key = TAIL_BASE + (received >> 1984)
    table = {r: table_word(r, *table_args) for r in [*keys, tail_key] if r >= 48}
    expected = table[tail_key] << 960
    for j, key in enumerate(keys):
        if key < 32:
            value = (received >> (64 * key)) & ((1 << 64) - 1)
        elif key < 48:
            value = 0
        else:
            value = table[key]
        lane = j % 4
        expected |= ((value >> (16 * lane)) & 65535) << (16 * j)
    machine = SampleMachine(bits(received, 2000), tables=table)
    assert replay(machine, program) == bits(expected, 1000)
    assert machine.ops == cost


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    bound, exponent = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
                       for c in manifest["claims"]]
    p, t = map(parity_table, (MAIN_COLUMNS, TAIL_COLUMNS))
    l, u = map(leaders, (MAIN_COLUMNS, TAIL_COLUMNS))
    a, b, h, k, stripe, prefix = success_counts(p, l, u)
    assert 0 < a < 1 and 0 < b < 1
    assert 1 - a ** 62 * b <= bound
    program = tuple(schedule())
    ops = Counter(op for op, *_ in program if op != "move")
    cost = sum(ops.values())
    assert cost == 283 and cost ** exponent.denominator <= 2 ** exponent.numerator
    rng = random.Random(2748)
    aliased = sum((j if j < 48 else TAIL_BASE + j) << (32 * j) for j in range(62))
    cases = [0, (1 << 2000) - 1, aliased, encode(rng.getrandbits(1000))]
    cases.extend(rng.getrandbits(2000) for _ in range(8))
    for y in cases:
        sample(y, program, (p, l, t, u), cost)
    print(json.dumps({"wide_histogram": h, "tail_histogram": k,
                      "stripe_loss_numerator": stripe, "prefix_loss_numerator": prefix,
                      "safe_wide_success": str(a), "tail_success": str(b),
                      "bound": str(bound), "ops": cost, "primitive_counts": ops,
                      "op_exponent": str(exponent), "sample_words": len(cases)}))


if __name__ == "__main__":
    main()
