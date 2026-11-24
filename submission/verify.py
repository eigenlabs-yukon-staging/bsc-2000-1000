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
    MAIN_COLUMNS, TAG_OFFSETS, TAIL_COLUMNS, TAIL_BASE, MAIN_BLOCKS, PLANE_GROUPS, encode, leaders, parity_table, schedule,
    table_word, quotient_choice,
)
from opmodel import Machine


def success_counts(main_parity, main_leaders, tail_leaders):
    histogram=Counter()
    tail_histogram=Counter(e.bit_count() for e in tail_leaders)
    multiplicity=Counter(main_parity)
    stripe=prefix=ordinary=0
    weights=[9**(32-d) for d in range(33)]
    for b in range(4):
        for s in range(1<<16):
            d,_=quotient_choice(s,b,main_leaders)
            histogram[d]+=1
            ordinary+=weights[d]
            stripe+=multiplicity[s^65535]*weights[d]
    for r in range(48):
        d,_=quotient_choice(main_parity[r],0,main_leaders)
        prefix+=weights[d]
    total=(1<<14)*ordinary-stripe-prefix
    safe_main=Fraction(total,(1<<16)*10**32)
    tail=Fraction(sum(c*9**(16-w) for w,c in tail_histogram.items()),10**16)
    return safe_main,tail,histogram,tail_histogram,stripe,prefix


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
    expected = 0
    for plane in range(PLANE_GROUPS):
        for lane in range(64):
            word = sum(((received >> (448 * plane + 64 * i + lane)) & 1) << i
                       for i in range(7))
            u = hamming_decode(word)
            for i in range(4):
                expected |= ((u >> i) & 1) << (256 * plane + 64 * i + lane)
    word_bits=[(expected>>(64*j))&((1<<64)-1) for j in range(16)]
    for j,key in enumerate(keys):
        if key<32:
            value=(received>>(64*key))&((1<<64)-1)
        elif key<48:
            value=0
        else:
            value=table[key]
        word,lane=divmod(j,4)
        word_bits[4*PLANE_GROUPS+word]|=(value<<(16*lane))&((1<<64)-1)
    word_bits[15]|=table[tail_key]
    word_bits[15]&=(1<<40)-1
    expected=sum(value<<(64*j) for j,value in enumerate(word_bits))
    machine = SampleMachine(bits(received, 2000), tables=table)
    assert replay(machine, program) == bits(expected, 1000)
    assert machine.ops == cost


def main():
    json.load(sys.stdin)
    manifest = json.loads((HERE / "submission.json").read_text())
    bound, exponent = [Fraction(int(c["score"]["numerator"]), int(c["score"]["denominator"]))
                       for c in manifest["claims"]]
    p, t = map(parity_table, (MAIN_COLUMNS, TAIL_COLUMNS))
    assert TAG_OFFSETS[0]==0 and TAG_OFFSETS[1]^TAG_OFFSETS[2]^TAG_OFFSETS[3]!=0
    l, u = leaders(MAIN_COLUMNS,16),leaders(TAIL_COLUMNS)
    a, b, h, k, stripe, prefix = success_counts(p, l, u)
    assert 0 < a < 1 and 0 < b < 1
    plane_success = hamming_success()
    assert 1 - plane_success ** (64 * PLANE_GROUPS) * a ** MAIN_BLOCKS * b <= bound
    program = tuple(schedule())
    ops = Counter(op for op, *_ in program if op != "move")
    cost = sum(ops.values())
    assert cost == 188 and cost ** exponent.denominator <= 2 ** exponent.numerator
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
