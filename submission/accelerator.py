"""Fixed wide-window field cube in scalar bsc-ops-v1 instructions."""
from bisect import bisect_right
from collections import Counter
from itertools import accumulate
from math import comb
import random

from construction import (
    ARRAYS, BUDGETS, GROUPS, HIGH_BITS, LIMIT, M, MASK, MASK64, MEMORY_END,
    MOD_LOW, RADIUS, WIDTHS, WORDS, carryless, cube, square,
)


SQUARE_WIDTH = 23 if M == 250 else 21
WINDOW_GROUPS = 2 if M == 250 else 3
WINDOW_WIDTH = 4 * WINDOW_GROUPS
REDUCTION_WORDS = int(M == 500)
SQUARE_WIDTHS = tuple(min(SQUARE_WIDTH, M - i) for i in range(0, M, SQUARE_WIDTH))
ROW_WORK = MEMORY_END
REDUCE_WIDE = ROW_WORK + WINDOW_GROUPS * 16 * WORDS
SQUARE_WIDE = REDUCE_WIDE + (1 << WINDOW_WIDTH) * REDUCTION_WORDS
SQUARE_BOUNDS = tuple(accumulate(
    ((1 << width) * WORDS for width in SQUARE_WIDTHS), initial=SQUARE_WIDE,
))
FAST_END = SQUARE_BOUNDS[-1]


def extra_table_word(address):
    if REDUCE_WIDE <= address < SQUARE_WIDE:
        index, word = divmod(address - REDUCE_WIDE, REDUCTION_WORDS)
        value = carryless(index, MOD_LOW)
    elif SQUARE_WIDE <= address < FAST_END:
        slot = bisect_right(SQUARE_BOUNDS, address) - 1
        index, word = divmod(address - SQUARE_BOUNDS[slot], WORDS)
        value = square(index << (SQUARE_WIDTH * slot))
    else:
        raise ValueError("outside fast fixed tables")
    return (value >> (64 * word)) & MASK64


def shift_reduce(regs, width):
    """Multiply a register vector by X**width modulo the field polynomial."""
    yield "shr", "c", regs[-1], HIGH_BITS - width
    for j in range(WORDS - 1, 0, -1):
        yield "shr", "k", regs[j - 1], 64 - width
        yield "shl", regs[j], regs[j], width
        yield "or_", regs[j], regs[j], "k"
    yield "shl", regs[0], regs[0], width
    yield "and_", regs[-1], regs[-1], (1 << HIGH_BITS) - 1
    if width == 1:
        yield "mul", "t", "c", MOD_LOW
        yield "xor", regs[0], regs[0], "t"
    elif M == 250:
        yield "mul", "t", "c", MOD_LOW
        yield "xor", regs[0], regs[0], "t"
        yield "shr", "t", "c", 5
        yield "xor", regs[1], regs[1], "t"
    else:
        yield "add", "a", "c", REDUCE_WIDE
        yield "load", "v", "a"
        yield "xor", regs[0], regs[0], "v"


def extract(regs, start, width):
    word, offset = divmod(start, 64)
    if offset:
        yield "shr", "n", regs[word], offset
    else:
        yield "move", "n", regs[word]
    if offset + width > 64:
        yield "shl", "t", regs[word + 1], 64 - offset
        yield "or_", "n", "n", "t"
    if offset + width != 64 and not (offset + width <= 64 and start + width == M):
        yield "and_", "n", "n", (1 << width) - 1


def accumulate_lookup(registers, row_base, initialize=False):
    yield "mul", "a", "n", len(registers)
    yield "add", "a", "a", row_base
    for j, reg in enumerate(registers):
        if j:
            yield "add", "t", "a", j
        yield "load", "v", "t" if j else "a"
        if initialize:
            yield "move", reg, "v"
        else:
            yield "xor", reg, reg, "v"


def cube_schedule(source, dest):
    """Compose a literal cube routine using six vectors and six scratch words."""
    x = tuple(f"x{j}" for j in range(WORDS))
    running = tuple(f"g{j}" for j in range(WORDS))
    powers = [tuple(f"b{i}_{j}" for j in range(WORDS)) for i in range(4)]

    for j, reg in enumerate(x):
        yield "load", reg, source + j
    for slot, width in enumerate(SQUARE_WIDTHS):
        yield from extract(x, SQUARE_WIDTH * slot, width)
        yield from accumulate_lookup(powers[0], SQUARE_BOUNDS[slot], slot == 0)

    previous = powers[0]
    for group in range(WINDOW_GROUPS):
        for i, registers in enumerate(powers):
            if group == i == 0:
                continue
            for dst, src in zip(registers, previous, strict=True):
                yield "move", dst, src
            yield from shift_reduce(registers, 1)
            previous = registers
        base = ROW_WORK + group * 16 * WORDS
        last = 0
        for step in range(1, 16):
            gray = step ^ (step >> 1)
            selected = powers[(gray ^ last).bit_length() - 1]
            for j, (reg, toggle) in enumerate(zip(running, selected, strict=True)):
                if step == 1:
                    yield "move", reg, toggle
                else:
                    yield "xor", reg, reg, toggle
                yield "store", None, base + gray * WORDS + j, reg
            last = gray

    windows = tuple(range(0, M, WINDOW_WIDTH))
    for start in reversed(windows):
        if start != windows[-1]:
            yield from shift_reduce(running, WINDOW_WIDTH)
        for i in range(WINDOW_GROUPS):
            position = start + 4 * i
            if position < M:
                yield from extract(x, position, min(4, M - position))
                yield from accumulate_lookup(running, ROW_WORK + 16 * WORDS * i,
                                             start == windows[-1] and i == 0)
    for j, reg in enumerate(running):
        yield "store", None, dest + j, reg


def operation_counts():
    return Counter(op for op, *_ in cube_schedule(ARRAYS["X0"], ARRAYS["C0"])
                   if op != "move")


def accelerated_work(caps):
    """All-input upper bound with the predecessor's fixed admission tables."""
    from dynamics import SCHEDULE, statistics
    bound = [0] * (LIMIT + 1)
    costs = statistics(10)
    ops = sum(operation_counts().values())
    for index in range(len(SCHEDULE) - 1, -1, -1):
        width = SCHEDULE[index]
        reserve = len(SCHEDULE) - 1 - index
        if GROUPS == 4:
            setup = 67 if width == 24 else 47
            step = 80 if width == 24 else 28
            if reserve == 0:
                step += ops + 200
            elif (index + 1) % len(WIDTHS) == 0:
                step += ops + 180
        else:
            setup = 90 if width == 24 else 47
            step = 100 if width == 24 else 28
            if reserve == 0:
                step += ops + 320
        continuation = [w + step for w in bound]
        worst = [0] * (LIMIT + 1)
        for (counts, _), effective in zip(costs[width], caps[index], strict=True):
            usage = [0] * (LIMIT + 1)
            for price, count in counts:
                for cap in range(price + reserve, LIMIT + 1):
                    usage[cap] += count * continuation[cap - price]
            for r, cap in enumerate(effective):
                worst[r] = max(worst[r], usage[cap])
        bound = [setup + w for w in worst]
    setup = 250 if GROUPS == 4 else 560
    assert len(BUDGETS) == RADIUS + 1
    return 1000 + sum(comb(M, d) * (44 + setup + ops + bound[cap])
                      for d, cap in enumerate(BUDGETS))


def counter_samples():
    from construction import table_word
    from opmodel import Machine
    rng = random.Random(611405)
    cases = [0, 1, MASK, 1 << (M - 1), *(rng.getrandbits(M) for _ in range(4))]
    source, dest = ARRAYS["X0"], ARRAYS["C0"]
    program = tuple(cube_schedule(source, dest))
    expected_ops = sum(operation_counts().values())
    assert len({dest for op, dest, *_ in program if op != "store"}) + 4 <= 64
    assert MEMORY_END == ROW_WORK < REDUCE_WIDE <= SQUARE_WIDE < FAST_END <= 1 << 32
    for value in cases:
        machine = Machine([0] * 2000)
        for j in range(WORDS):
            machine.mem[source + j] = (value >> (64 * j)) & MASK64
        regs = {}
        for op, register, *arguments in program:
            args = [regs[a] if isinstance(a, str) else a for a in arguments]
            if op == "move":
                regs[register] = args[0]
            elif op == "store":
                machine.store(*args)
            else:
                if op == "load" and REDUCE_WIDE <= args[0] < FAST_END:
                    if args[0] not in machine.mem:
                        machine.mem[args[0]] = table_word(args[0])
                regs[register] = getattr(machine, op)(*args)
        result = sum(machine.mem[dest + j] << (64 * j) for j in range(WORDS))
        assert result == cube(value)
        assert machine.ops == expected_ops
    return {"cube_ops": expected_ops, "cube_replays": len(cases),
            "cube_counts": operation_counts(), "fast_memory_end": FAST_END}
