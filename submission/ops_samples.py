"""Replay literal decoder schedules using only exact sparse table snapshots."""
from collections import Counter

from construction import OFFSET, SYNDROME_COLUMNS, TAIL_BASE, encode, schedule, table_word, tail_leaders
from opmodel import Machine


def bits(word, width):
    return [(word >> i) & 1 for i in range(width)]


def replay(y, expected, program, keys, leaders):
    tails = tail_leaders()
    tables = {key: table_word(key, leaders, tails) for key in keys}
    machine = Machine(bits(y, 2000), tables=tables)
    regs = {}
    for op, dest, *params in program:
        args = [regs[a] if isinstance(a, str) else a for a in params]
        if op == "store":
            machine.store(*args)
        else:
            if op == "load" and args[0] >= 48:
                assert args[0] in tables
            regs[dest] = getattr(machine, op)(*args)
    assert machine.output_bits() == bits(expected, 1000)
    assert machine.ops == 283


def check_runtime(exponent):
    program = tuple(schedule())
    counts = Counter(op for op, *_ in program)
    assert sum(counts.values()) == 283
    assert 283**exponent.denominator <= 2**exponent.numerator
    replay(0, 0, program, {TAIL_BASE}, {})
    reserved = sum(47 << (64 * j + 32) for j in range(1, 31))
    replay(reserved, 0, program, {TAIL_BASE}, {})

    expected = sum(OFFSET << (32 * j) for j in range(31)) | (195 << 992)
    y = encode(expected)
    leaders = {0: 0}
    keys = {OFFSET, TAIL_BASE + (y >> 1984)}
    for j in range(31):
        bit = j % 16
        s = 1 << bit
        key = OFFSET ^ s
        assert 48 <= key < TAIL_BASE and 48 <= OFFSET < TAIL_BASE
        keys.add(key)
        assert SYNDROME_COLUMNS.count(s) == 1
        leaders[s] = 1 << (32 + bit)
        y ^= 1 << (64 * j + 32 + bit)
    replay(y, expected, program, keys, leaders)
    return {"ops": 283, "counts": counts, "exact_sample_replays": 3}
