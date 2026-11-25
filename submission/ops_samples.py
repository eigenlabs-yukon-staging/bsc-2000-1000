"""Replay four exact sparse snapshots on the unchanged organizer counter."""
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
    assert machine.ops == 279


def check_runtime(exponent):
    program = tuple(schedule())
    counts = Counter(op for op, *_ in program)
    assert sum(counts.values()) == 279
    assert 279**exponent.denominator <= 2**exponent.numerator

    # In the reserved snapshot all first reads use low half zero.
    # Input word zero cancels the marked wide keys to zero, while
    # the tail's shifted marked key is in initial output word 47.
    reserved = sum(47 << (64*j + 32) for j in range(30))
    reserved |= 47 << (64*30 + 40)
    run_case(0, 0, program, ())
    run_case(reserved, 0, program, ())

    payload = sum(OFFSET << (32*i) for i in range(31)) | (195 << 992)
    received = encode(payload)
    table_addresses = {OFFSET}
    for j in range(30):
        received ^= 1 << (64*j + 32 + j)
        address = OFFSET ^ (1 << j)
        assert address >= TAIL_LIMIT
        table_addresses.add(address)
    received ^= 1 << (64*30 + 40)
    short_address = (OFFSET >> 8) ^ 1
    assert 48 <= short_address < TAIL_LIMIT
    table_addresses.add(short_address)
    assert ((payload >> 960) & MASK32) == OFFSET
    run_case(received, payload, program, table_addresses)

    # This low bank address gets full wide parity and wide correction.
    # Exercise it first as a parity key and then as a correction key.
    restored = 0x85c99d
    wide_noise = 0x1000000020000002
    mixed_payload = (payload & ~MASK32) | restored
    mixed_received = encode(mixed_payload) ^ (wide_noise << 64)
    mixed_addresses = {restored, OFFSET, OFFSET ^ (wide_noise & MASK32), OFFSET >> 8}
    run_case(mixed_received, mixed_payload, program, mixed_addresses)

    return {"ops": 279, "counts": counts, "exact_sample_replays": 4}
