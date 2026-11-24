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
    MAIN_COLUMNS, TAIL_COLUMNS, TAIL_BASE, MAIN_BLOCKS, PLANE_GROUPS, encode, leaders, parity_table, schedule,
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


def hamming_success():
    cols = (3, 5, 6, 7, 1, 2, 4)
    leaders = {0: 0} | {c: 1 << i for i, c in enumerate(cols)}
    assert len(leaders) == 8
    return Fraction(sum(9 ** (7 - e.bit_count()) for e in leaders.values()), 10 ** 7)


def hamming_decode(received):
    cols = (3, 5, 6, 7, 1, 2, 4)
    syndrome = 0
    for i, col in enumerate(cols):
        if (received >> i) & 1:
            syndrome ^= col
    flip = (1 << cols.index(syndrome)) if syndrome else 0
    data = (received ^ flip) & 15
    return data ^ (7 if data & 8 else 0)


def sample(received, program, table_args, cost):
    keys = [(received >> (448 * PLANE_GROUPS + 32 * j)) & 0xFFFFFFFF for j in range(MAIN_BLOCKS)]
    tail_key = TAIL_BASE + (received >> 1984)
    table = {r: table_word(r, *table_args) for r in [*keys, tail_key] if r >= 48}
    expected = table[tail_key] << 960
    for plane in range(PLANE_GROUPS):
        for lane in range(64):
            word = sum(((received >> (448 * plane + 64 * i + lane)) & 1) << i
                       for i in range(7))
            u = hamming_decode(word)
            for i in range(4):
                expected |= ((u >> i) & 1) << (256 * plane + 64 * i + lane)
    for j, key in enumerate(keys):
        if key < 32:
            value = (received >> (64 * key)) & ((1 << 64) - 1)
        elif key < 48:
            value = 0
        else:
            value = table[key]
        lane = j % 4
        expected |= ((value >> (16 * lane)) & 65535) << (256 * PLANE_GROUPS + 16 * j)
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
    plane_success = hamming_success()
    assert 1 - plane_success ** (64 * PLANE_GROUPS) * a ** MAIN_BLOCKS * b <= bound
    program = tuple(schedule())
    ops = Counter(op for op, *_ in program if op != "move")
    cost = sum(ops.values())
    assert cost == 195 and cost ** exponent.denominator <= 2 ** exponent.numerator
    rng = random.Random(2748)
    special_keys = [0, 1, 31, 32, 33, 47, 48]
    special_keys += [TAIL_BASE + j for j in range(MAIN_BLOCKS - len(special_keys))]
    aliased = sum(key << (448 * PLANE_GROUPS + 32 * j)
                  for j, key in enumerate(special_keys))
    cases = [0, (1 << 2000) - 1, aliased, encode(rng.getrandbits(1000))]
    cases.extend(rng.getrandbits(2000) for _ in range(8))
    for y in cases:
        sample(y, program, (p, l, t, u), cost)
    print(json.dumps({"plane_success": str(plane_success), "wide_histogram": h, "tail_histogram": k,
                      "stripe_loss_numerator": stripe, "prefix_loss_numerator": prefix,
                      "safe_wide_success": str(a), "tail_success": str(b),
                      "bound": str(bound), "ops": cost, "primitive_counts": ops,
                      "op_exponent": str(exponent), "sample_words": len(cases)}))


if __name__ == "__main__":
    main()
