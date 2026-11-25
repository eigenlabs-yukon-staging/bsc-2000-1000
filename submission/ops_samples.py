"""Replay three exact sparse snapshots on the unchanged organizer counter."""
from collections import Counter

from construction import MASK32, OFFSET, TAIL_LIMIT, encode, schedule, table_word
from opmodel import Machine


def bits(word, width):
    return [(word >> i) & 1 for i in range(width)]


def run_case(received, expected, program, possible_table_addresses):
    # All nonreserved addresses touched on these snapshots are precomputed.
    tables = {address: table_word(address) for address in possible_table_addresses}
    machine = Machine(bits(received, 2000), tables=tables)
    registers = {}
    for op, dest, *args in program:
        values = [registers[arg] if isinstance(arg, str) else arg for arg in args]
        if op == "store":
            machine.store(*values)
        else:
            if op == "load" and values[0] >= 48:
                assert values[0] in tables
            registers[dest] = getattr(machine, op)(*values)
    assert machine.output_bits() == bits(expected, 1000)
    assert machine.ops == 147


def check_runtime(exponent):
    program = tuple(schedule())
    counts = Counter(op for op, *_ in program)
    assert sum(counts.values()) == 147
    assert 147**exponent.denominator <= 2**exponent.numerator

    # Input 0 is zero and every wide/shortened low part is zero. Their
    # parity-only changes send their second reads into output word 47,
    # which is still zero while all lanes are assembled.
    reserved = sum(47 << (64 * j + 32) for j in range(21, 27))
    reserved |= 47 << (64 * 27 + 40)
    run_case(0, 0, program, ())
    run_case(reserved, 0, program, ())

    low768 = int.from_bytes(bytes(range(96)), "little")
    payload = low768 | sum(OFFSET << (768 + 32 * i) for i in range(7))
    payload |= 195 << 992
    received = encode(payload)
    # A single noisy Hamming plane affects all its lanes by one bit each.
    received ^= ((1 << 64) - 1) << (7 * 64)
    table_addresses = {OFFSET}
    for j in range(6):
        received ^= 1 << (64 * (21 + j) + 32 + j)
        addr = OFFSET ^ (1 << j)
        assert addr >= TAIL_LIMIT
        table_addresses.add(addr)
    received ^= 1 << (64 * 27 + 40)
    short_addr = (OFFSET >> 8) ^ 1
    assert 48 <= short_addr < TAIL_LIMIT
    table_addresses.add(short_addr)
    assert ((payload >> 960) & MASK32) == OFFSET
    run_case(received, payload, program, table_addresses)

    return {"ops": 147, "counts": counts, "exact_sample_replays": 3}
