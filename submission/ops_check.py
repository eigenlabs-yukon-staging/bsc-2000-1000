"""Validate the primitive list, table addressing, and packing with the reference counter."""
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path
import random
import sys

HERE = Path("/snapshot/files") if Path("/snapshot/files/construction.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from construction import (
    LAYOUT, block_columns, cap, encode, parity, schedule, syndrome_columns,
    table_addresses,
)
from opmodel import Machine


def low_weight_leaders(k):
    cols = syndrome_columns(k)
    leaders = {}
    for w in range(3):
        for positions in itertools.combinations(range(2 * k), w):
            syndrome = 0
            for j in positions:
                syndrome ^= cols[j]
            leaders.setdefault(syndrome, sum(1 << j for j in positions))
    return leaders


def sampled_preload(received, leaders):
    """Compute exact preloaded values for the addresses used in this sample."""
    offset = 0
    addresses, _ = table_addresses()
    tables = {}
    result = 0
    for k in LAYOUT:
        p_base, c_base = addresses[k]
        a = (received >> (2 * offset)) & cap(k)
        v = (received >> (2 * offset + k)) & cap(k)
        p_a = parity(a, block_columns(k))
        c = leaders[k][p_a ^ v] & cap(k)
        tables[p_base + a] = p_a
        tables[c_base + (p_a ^ v)] = c
        result |= (a ^ c) << offset
        offset += k
    return tables, result


class CheckedMachine(Machine):
    def load(self, address):
        if address >= 48 and address not in self.mem:
            raise ValueError("missing preloaded table entry for a sample")
        return super().load(address)


def replay(machine, program):
    registers = {}
    for op, dst, *args in program:
        values = [registers[a] if isinstance(a, str) else a for a in args]
        if op == "move":
            registers[dst] = values[0]
        elif op == "store":
            machine.store(*values)
        else:
            registers[dst] = getattr(machine, op)(*values)
    return machine.output_bits()


def bits(value, width):
    return [(value >> i) & 1 for i in range(width)]


def main():
    json.load(sys.stdin)
    score = json.loads((HERE / "submission.json").read_text())["claims"][1]["score"]
    exponent = Fraction(int(score["numerator"]), int(score["denominator"]))
    program = list(schedule())
    counts = Counter(op for op, *_ in program if op != "move")
    cost = sum(counts.values())
    assert cost == 519
    assert cost ** exponent.denominator <= 2 ** exponent.numerator
    _, end = table_addresses()
    assert end <= 1 << 32
    rng = random.Random(419)
    leaders = {k: low_weight_leaders(k) for k in (30, 10)}
    cases = 0
    for i in range(20):
        message = (0 if i == 0 else (1 << 1000) - 1 if i == 1 else rng.getrandbits(1000))
        received = encode(message)
        if i >= 2:
            offset = 0
            for k in LAYOUT:
                for j in rng.sample(range(2 * k), i % 3):
                    received ^= 1 << (offset + j)
                offset += 2 * k
        tables, expected = sampled_preload(received, leaders)
        machine = CheckedMachine(bits(received, 2000), tables=tables)
        actual = replay(machine, program)
        assert machine.ops == cost and actual == bits(expected, 1000)
        if i < 2:
            assert actual == bits(message, 1000)
        cases += 1
    print(json.dumps({"op_counts": counts, "ops": cost, "exponent": str(exponent),
                      "integer_power_comparison": True, "memory_words": end,
                      "sample_words": cases}))


if __name__ == "__main__":
    main()
