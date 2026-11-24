"""Verify the extended Golay construction and its fixed RAM program."""
import itertools
import json
from pathlib import Path
import random
import sys
from array import array
from collections import Counter
from fractions import Fraction

HERE = Path("/snapshot/files") if Path("/snapshot/files/opmodel.py").is_file() else Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from opmodel import Machine

GENERATOR = sum(1 << i for i in (0, 2, 4, 5, 6, 10, 11))
BASE = 48
MASK = (1 << 24) - 1


def codeword(message):
    body = 0
    for i in range(12):
        if message >> i & 1:
            body ^= GENERATOR << i
    return body | ((body.bit_count() & 1) << 23)


def syndrome(word):
    parity = word.bit_count() & 1
    body = word & ((1 << 23) - 1)
    for i in range(22, 10, -1):
        if body >> i & 1:
            body ^= GENERATOR << (i - 11)
    return body | (parity << 11)


def leaders():
    result = {}
    for weight in range(5):
        for positions in itertools.combinations(range(24), weight):
            noise = sum(1 << i for i in positions)
            result.setdefault(syndrome(noise), noise)
    return result


def table(chosen):
    result = array("H", [65535]) * (1 << 24)
    for message in range(1 << 12):
        sent = codeword(message)
        for noise in chosen.values():
            received = sent ^ noise
            assert result[received] == 65535
            result[received] = message
    return result


def schedule():
    for j in range(32):
        yield "load", f"Y{j}", j
    for j in range(16):
        yield "move", f"U{j}", 0
    for i in range(83):
        word, offset = divmod(24 * i, 64)
        yield ("shr", "T", f"Y{word}", offset) if offset else ("move", "T", f"Y{word}")
        if offset > 40:
            yield "shl", "S", f"Y{word+1}", 64-offset
            yield "or_", "T", "T", "S"
        if offset != 40:
            yield "and_", "T", "T", MASK
        yield "add", "T", "T", BASE
        yield "load", "D", "T"
        word, offset = divmod(12 * i, 64)
        if offset:
            yield "shl", "S", "D", offset
            yield "or_", f"U{word}", f"U{word}", "S"
        else:
            yield "move", f"U{word}", "D"
        if offset > 52:
            yield "shr", f"U{word+1}", "D", 64-offset
    yield "shr", "T", "Y31", 8
    yield "and_", "T", "T", 15
    yield "shl", "T", "T", 36
    yield "or_", "U15", "U15", "T"
    for j in range(16):
        yield "store", None, 32+j, f"U{j}"


class TableMemory(dict):
    def __init__(self, existing, lookup):
        super().__init__(existing)
        self.lookup = lookup

    def get(self, address, default=None):
        if BASE <= address < BASE + len(self.lookup):
            return self.lookup[address - BASE]
        return super().get(address, default)


def execute(machine, program):
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


def encode(message):
    encoded = sum(codeword((message >> (12*i)) & 4095) << (24*i) for i in range(83))
    tail = message >> 996
    return encoded | (tail << 1992) | (tail << 1996)


def decode_reference(received, chosen):
    result = 0
    for i in range(83):
        block = (received >> (24*i)) & MASK
        corrected = block ^ chosen[syndrome(block)]
        for j in range(22, 10, -1):
            if corrected >> j & 1:
                result |= 1 << (12*i+j-11)
                corrected ^= GENERATOR << (j-11)
    return result | (((received >> 1992) & 15) << 996)


def bits(value, size):
    return [(value >> i) & 1 for i in range(size)]


def exact_values(chosen, program):
    counts = Counter(noise.bit_count() for noise in chosen.values())
    assert counts == {0: 1, 1: 24, 2: 276, 3: 2024, 4: 1771}
    component = sum(count * 9**(24-weight) for weight, count in counts.items())
    error = 1 - Fraction(component**83 * 9**4, 10**1996)
    cost = sum(op != "move" for op, *_ in program)
    return counts, component, error, cost


def main():
    json.load(sys.stdin)
    metadata = json.loads((HERE / "submission.json").read_text())
    bound, exponent = [
        Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
        for c in metadata["claims"]
    ]
    chosen = leaders()
    program = list(schedule())
    counts, component, error, cost = exact_values(chosen, program)
    assert error <= bound
    assert cost**exponent.denominator <= 2**exponent.numerator
    lookup = table(chosen)
    rng = random.Random(173)
    messages = [0, (1 << 1000) - 1] + [rng.getrandbits(1000) for _ in range(6)]
    words = [encode(message) for message in messages]
    words.extend([0, (1 << 2000) - 1] + [rng.getrandbits(2000) for _ in range(6)])
    for received in words:
        machine = Machine(bits(received, 2000))
        machine.mem = TableMemory(machine.mem, lookup)
        result = execute(machine, program)
        assert machine.ops == cost
        assert result == bits(decode_reference(received, chosen), 1000)
    assert [decode_reference(w, chosen) for w in words[:len(messages)]] == messages
    print(json.dumps({"leader_histogram": counts, "component_success_numerator": component,
                      "bound": str(bound), "bound_verified": True, "ops": cost,
                      "exponent": str(exponent), "integer_power_comparison": True,
                      "table_entries": len(lookup), "sample_words": len(words)}))


if __name__ == "__main__":
    main()
