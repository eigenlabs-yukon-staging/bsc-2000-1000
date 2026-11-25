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
    assert machine.ops == 191


def check_runtime(exponent):
    program = tuple(schedule())
    counts = Counter(op for op, *_ in program)
    assert sum(counts.values()) == 191
    assert 191**exponent.denominator <= 2**exponent.numerator

    # Word zero is zero, so the exceptional wide and tail syndrome
    # keys point to initial output word 47 before the final stores.
    reserved = sum(47 << (64*j+32) for j in range(14,28))
    reserved |= 47 << (64*28+40)
    run_case(0, 0, program, ())
    run_case(reserved, 0, program, ())

    low512 = int.from_bytes(bytes(range(64)), "little")
    payload = low512 | sum(OFFSET << (512+32*j) for j in range(15)) | (195 << 992)
    received = encode(payload)
    for group in range(2):
        received ^= ((1 << 64)-1) << (64*(7*group+3))
    addresses = {OFFSET}
    for j in range(14):
        received ^= 1 << (64*(14+j)+32+j)
        address = OFFSET ^ (1 << j)
        assert address >= TAIL_LIMIT
        addresses.add(address)
    received ^= 1 << (64*28+40)
    addresses.add((OFFSET >> 8) ^ 1)
    run_case(received, payload, program, addresses)

    restored = 0x85c99d
    wide_noise = 0x1000000020000002
    mixed_payload = (payload & ~(MASK32 << 512)) | (restored << 512)
    mixed_received = encode(mixed_payload) ^ (wide_noise << (64*15))
    mixed_addresses = {restored, OFFSET, OFFSET ^ (wide_noise & MASK32), OFFSET >> 8}
    run_case(mixed_received, mixed_payload, program, mixed_addresses)

    return {"ops": 191, "counts": counts, "exact_sample_replays": 4}
