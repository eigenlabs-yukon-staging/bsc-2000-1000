"""Reuse the first field cube and its derivative within a fixed search prefix."""
from bisect import bisect_right
from functools import lru_cache
from itertools import accumulate
from math import comb
import random

from accelerator import FAST_END, accumulate_lookup, extract, operation_counts, shift_reduce
from construction import (
    ANCHOR_BITS, ARRAYS, BUDGETS, LIMIT, M, MASK64, OFFSETS, WORDS,
    cube, remainder, square, table_word,
)

FIRST_TAIL_CELL = 15
TAIL_START = OFFSETS[FIRST_TAIL_CELL]
TAIL_WIDTH = M - TAIL_START
TAIL_WORD = TAIL_START // 64
TAIL_OFFSET = TAIL_START % 64

MODE_STATE = FAST_END
READY_STATE = MODE_STATE + 1
BASE_CUBE = READY_STATE + 1
BASE_TAIL = BASE_CUBE + WORDS
CROSS_WORK = BASE_TAIL + WORDS - TAIL_WORD


def widths(bits, width):
    return tuple(min(width, bits - j) for j in range(0, bits, width))


def offsets(parts):
    return tuple(accumulate(parts, initial=0))


def bases(start, row_sizes):
    return tuple(accumulate(row_sizes, initial=start))


CROSS_WIDTHS = widths(TAIL_WIDTH, 3)
CROSS_OFFSETS = offsets(CROSS_WIDTHS)
CROSS_BASES = bases(CROSS_WORK, ((1 << w) * WORDS for w in CROSS_WIDTHS))
MODE_BASE = CROSS_BASES[-1]
SQUARE_BYTE_BASE = MODE_BASE + LIMIT + 1
LINEAR_BASE = SQUARE_BYTE_BASE + 4 * 256
LINEAR_WIDTHS = widths(M, 17)
LINEAR_OFFSETS = offsets(LINEAR_WIDTHS)
LINEAR_BASES = bases(LINEAR_BASE, ((1 << w) * 2 * WORDS for w in LINEAR_WIDTHS))
RAW_BITS = 3 * TAIL_WIDTH - 2
RAW_WIDTHS = widths(RAW_BITS, 19)
RAW_OFFSETS = offsets(RAW_WIDTHS)
RAW_BASES = bases(LINEAR_BASES[-1], ((1 << w) * WORDS for w in RAW_WIDTHS))
PREFIX_END = RAW_BASES[-1]
RESIDUE32 = tuple(0x11111111 << i for i in range(4))
RESIDUE64 = tuple(0x1111111111111111 << i for i in range(4))

def prefix_word(address):
    if MODE_BASE <= address < SQUARE_BYTE_BASE:
        return work_and_modes()[1][address - MODE_BASE]
    if SQUARE_BYTE_BASE <= address < LINEAR_BASE:
        slot, byte = divmod(address - SQUARE_BYTE_BASE, 256)
        return sum(((byte >> j) & 1) << (16 * slot + 2 * j) for j in range(8))
    if LINEAR_BASE <= address < LINEAR_BASES[-1]:
        j = bisect_right(LINEAR_BASES, address) - 1
        key, word = divmod(address - LINEAR_BASES[j], 2 * WORDS)
        part = key << LINEAR_OFFSETS[j]
        if word < WORDS:
            value = remainder(square(part) << TAIL_START)
        else:
            word -= WORDS
            value = remainder(part << (2 * TAIL_START))
    elif RAW_BASES[0] <= address < PREFIX_END:
        j = bisect_right(RAW_BASES, address) - 1
        key, word = divmod(address - RAW_BASES[j], WORDS)
        value = remainder(key << (3 * TAIL_START + RAW_OFFSETS[j]))
    else:
        raise ValueError("outside fixed prefix tables")
    return (value >> (64 * word)) & MASK64


def cache_setup(source, dest):
    """Save the first cube and setup per-prefix positional derivative rows."""
    for j in range(WORDS):
        yield "load", "v", dest + j
        yield "store", None, BASE_CUBE + j, "v"
    for j in range(TAIL_WORD, WORDS):
        yield "load", "v", source + j
        yield "store", None, BASE_TAIL + j - TAIL_WORD, "v"

    running = tuple(f"r{j}" for j in range(WORDS))
    v = tuple(f"v{j}" for j in range(WORDS))
    w = tuple(f"w{j}" for j in range(WORDS))
    basis = [tuple(f"q{i}_{j}" for j in range(WORDS)) for i in range(3)]
    for j, reg in enumerate(running):
        yield "load", reg, source + j
    for j, width in enumerate(LINEAR_WIDTHS):
        yield from extract(running, LINEAR_OFFSETS[j], width)
        yield from accumulate_lookup(v + w, LINEAR_BASES[j], j == 0)

    for group, width in enumerate(CROSS_WIDTHS):
        for local in range(width):
            for j in range(WORDS):
                yield "xor", basis[local][j], v[j], w[j]
            if CROSS_OFFSETS[group] + local + 1 < TAIL_WIDTH:
                yield from shift_reduce(v, 1)
                yield from shift_reduce(w, 2)
        last = 0
        for step in range(1, 1 << width):
            key = step ^ (step >> 1)
            chosen = basis[(key ^ last).bit_length() - 1]
            for j, reg in enumerate(running):
                if step == 1:
                    yield "move", reg, chosen[j]
                else:
                    yield "xor", reg, reg, chosen[j]
                yield "store", None, CROSS_BASES[group] + key * WORDS + j, reg
            last = key


def raw_cube(diff, raw):
    """Compute the ordinary binary-polynomial cube of the normalized tail."""
    a = [tuple(f"a{i}_{r}" for r in range(4)) for i in range((TAIL_WIDTH + 31) // 32)]
    b = [tuple(f"b{i}_{r}" for r in range(2)) for i in range((2 * TAIL_WIDTH + 30) // 32)]
    for i, terms in enumerate(a):
        width = min(32, TAIL_WIDTH - 32 * i)
        yield from extract(diff, TAIL_OFFSET + 32 * i, width)
        yield "move", "limb", "n"
        for r, reg in enumerate(terms):
            yield "and_", reg, "limb", RESIDUE32[r]
        for j in range((width + 7) // 8):
            if j:
                yield "shr", "n", "limb", 8 * j
            if 8 * (j + 1) < width:
                yield "and_", "n", "n", 255
            yield "add", "a", "n", SQUARE_BYTE_BASE + j * 256
            yield "load", "v", "a"
            if j:
                yield "or_", "sq", "sq", "v"
            else:
                yield "move", "sq", "v"
        for part in range(2):
            if 2 * i + part >= len(b):
                break
            if part:
                yield "shr", "sq", "sq", 32
            yield "and_", b[2 * i + part][0], "sq", RESIDUE32[0]
            yield "and_", b[2 * i + part][1], "sq", RESIDUE32[2]

    for reg in raw:
        yield "move", reg, 0
    for i, terms in enumerate(a):
        for j, (even, odd) in enumerate(b):
            for r in range(4):
                yield "mul", "p", terms[r], even
                yield "mul", "q", terms[(r - 2) % 4], odd
                yield "xor", "p", "p", "q"
                yield "and_", "p", "p", RESIDUE64[r]
                if r:
                    yield "xor", "prod", "prod", "p"
                else:
                    yield "move", "prod", "p"
            word, half = divmod(i + j, 2)
            if half:
                yield "shl", "p", "prod", 32
                yield "xor", raw[word], raw[word], "p"
                if word + 1 < len(raw):
                    yield "shr", "p", "prod", 32
                    yield "xor", raw[word + 1], raw[word + 1], "p"
            else:
                yield "xor", raw[word], raw[word], "prod"


def cache_reuse(source, dest):
    """Recover the current cube from the first cube and a tail-only difference."""
    diff = tuple(f"d{j}" for j in range(WORDS - TAIL_WORD))
    running = tuple(f"a{i}_{j}" for i in range(2) for j in range(4))
    raw = tuple(f"z{j}" for j in range((RAW_BITS + 63) // 64))
    for j, reg in enumerate(diff):
        yield "load", reg, source + TAIL_WORD + j
        yield "load", "v", BASE_TAIL + j
        yield "xor", reg, reg, "v"
    yield from raw_cube(diff, raw)
    for j, reg in enumerate(running):
        yield "load", reg, BASE_CUBE + j
    for j, width in enumerate(CROSS_WIDTHS):
        yield from extract(diff, TAIL_OFFSET + CROSS_OFFSETS[j], width)
        yield from accumulate_lookup(running, CROSS_BASES[j])
    for j, width in enumerate(RAW_WIDTHS):
        yield from extract(raw, RAW_OFFSETS[j], width)
        yield from accumulate_lookup(running, RAW_BASES[j])
    for j, reg in enumerate(running):
        yield "store", None, dest + j, reg


def count(schedule):
    return sum(op != "move" for op, *_ in schedule(ARRAYS["X1"], ARRAYS["C1"]))


def step_bound(bound, counts, effective, reserve, setup, step):
    continuation = [c + step for c in bound]
    worst = [0] * (LIMIT + 1)
    for (profile, _), row in zip(counts, effective, strict=True):
        usage = [0] * (LIMIT + 1)
        for price, n in profile:
            for cap in range(price + reserve, LIMIT + 1):
                usage[cap] += n * continuation[cap - price]
        for r, cap in enumerate(row):
            worst[r] = max(worst[r], usage[cap])
    return [setup + value for value in worst]


@lru_cache(None)
def work_and_modes():
    from dynamics import SCHEDULE, derive, statistics
    _, caps, _ = derive(10, True)
    counts = statistics(10)
    full_cube = sum(operation_counts().values())
    prep = full_cube + count(cache_setup)
    reuse = count(cache_reuse)
    full = saved = [0] * (LIMIT + 1)
    mode = None
    for j in range(len(SCHEDULE) - 1, -1, -1):
        width = SCHEDULE[j]
        reserve = len(SCHEDULE) - j - 1
        h, g = (90, 100) if width == 24 else (47, 28)
        full_step = g + (320 + full_cube + 3 if reserve == 0 else 0)
        reuse_step = g + (320 + reuse + 6 if reserve == 0 else 0)
        full = step_bound(full, counts[width], caps[j], reserve, h, full_step)
        if j >= FIRST_TAIL_CELL:
            saved = step_bound(saved, counts[width], caps[j], reserve, h, reuse_step)
            if j == FIRST_TAIL_CELL:
                mode = tuple(int(s + prep < a) for a, s in zip(full, saved))
                full = [4 + min(a, s + prep) for a, s in zip(full, saved)]
    result = 1000 + sum(comb(ANCHOR_BITS, d) * (44 + 560 + full_cube + full[cap])
                        for d, cap in enumerate(BUDGETS))
    return result, mode


def cache_samples():
    from opmodel import Machine
    assert M == 500 and TAIL_START == 360 and TAIL_WIDTH == 140
    assert FAST_END < MODE_BASE < LINEAR_BASE < PREFIX_END <= 1 << 32
    source, dest = ARRAYS["X1"], ARRAYS["C1"]
    schedules = (tuple(cache_setup(source, dest)), tuple(cache_reuse(source, dest)))
    for program in schedules:
        assert len({dst for op, dst, *_ in program if op != "store"}) + 4 <= 64

    def replay(machine, program):
        regs = {}
        for op, dst, *params in program:
            args = [regs[a] if isinstance(a, str) else a for a in params]
            if op == "move":
                regs[dst] = args[0]
            elif op == "store":
                machine.store(*args)
            else:
                if op == "load" and args[0] not in machine.mem:
                    if SQUARE_BYTE_BASE <= args[0] < PREFIX_END or args[0] < FAST_END and args[0] >= 48:
                        machine.mem[args[0]] = table_word(args[0])
                regs[dst] = getattr(machine, op)(*args)

    def install(machine, address, value):
        for j in range(WORDS):
            machine.mem[address+j] = (value >> (64*j)) & MASK64

    rng = random.Random(255812)
    bases = (0, (1 << M) - 1, *(rng.getrandbits(M) for _ in range(3)))
    deltas = (0, 1, (1 << TAIL_WIDTH) - 1, rng.getrandbits(TAIL_WIDTH))
    results = 0
    for base in bases:
        machine = Machine([0] * 2000)
        install(machine, source, base)
        install(machine, dest, cube(base))
        replay(machine, schedules[0])
        assert machine.ops == count(cache_setup)
        for d in deltas:
            v = base ^ (d << TAIL_START)
            install(machine, source, v)
            before = machine.ops
            replay(machine, schedules[1])
            observed = sum(machine.mem[dest+j] << (64*j) for j in range(WORDS))
            assert observed == cube(v)
            assert machine.ops - before == count(cache_reuse)
            results += 1
    w, mode = work_and_modes()
    assert len(mode) == LIMIT + 1
    return {"setup_ops": count(cache_setup), "reuse_ops": count(cache_reuse),
            "reuse_replays": results, "mode_count": sum(mode),
            "prefix_memory_end": PREFIX_END, "W": str(w)}
