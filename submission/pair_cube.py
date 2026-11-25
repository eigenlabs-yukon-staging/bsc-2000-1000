"""Fixed schedule for a degree-50 field cube with six pair tables."""
from bisect import bisect_right
from collections import Counter
from itertools import accumulate, combinations
from math import comb
import random

from construction import ARRAYS, BUDGETS, COMPONENTS, LIMIT, M, MASK, MEMORY_END, RADIUS, WIDTHS, cube


CHUNKS = (13, 13, 13, 11)
OFFSETS = tuple(accumulate(CHUNKS, initial=0))
PAIRS = tuple(combinations(range(len(CHUNKS)), 2))
UNARY_OWNER = (1, 0, 0, 0)
BASES = tuple(accumulate((1 << (CHUNKS[i] + CHUNKS[j]) for i, j in PAIRS),
                         initial=MEMORY_END))
PAIR_END = BASES[-1]


def pair_word(address):
    if not MEMORY_END <= address < PAIR_END:
        raise ValueError("outside the pair cube tables")
    index = bisect_right(BASES, address) - 1
    i, j = PAIRS[index]
    key = address - BASES[index]
    a = (key & ((1 << CHUNKS[i]) - 1)) << OFFSETS[i]
    b = (key >> CHUNKS[i]) << OFFSETS[j]
    value = cube(a ^ b)
    if UNARY_OWNER[i] != j:
        value ^= cube(a)
    if UNARY_OWNER[j] != i:
        value ^= cube(b)
    return value


def cube_schedule(source, dest):
    yield "load", "x", source
    chunks = tuple(f"part{i}" for i in range(len(CHUNKS)))
    raised = tuple(f"high{i}" for i in range(len(CHUNKS)))
    for i, (part, width, start) in enumerate(zip(chunks, CHUNKS, OFFSETS[:-1], strict=True)):
        if start:
            yield "shr", part, "x", start
        else:
            yield "move", part, "x"
        if start + width < M:
            yield "and_", part, part, (1 << width) - 1
        if i:
            yield "shl", raised[i], part, 13
    for idx, (i, j) in enumerate(PAIRS):
        yield "or_", "key", chunks[i], raised[j]
        yield "add", "address", "key", BASES[idx]
        yield "load", "value", "address"
        if idx == 0:
            yield "move", "result", "value"
        else:
            yield "xor", "result", "result", "value"
    yield "store", None, dest, "result"


def operation_counts():
    return Counter(op for op, *_ in cube_schedule(ARRAYS["X0"], ARRAYS["C0"])
                   if op != "move")


def state_work(caps):
    from dynamics import SCHEDULE, statistics

    continuation = [0] * (LIMIT + 1)
    rows = statistics(10)
    cube_ops = sum(operation_counts().values())
    for i in range(len(SCHEDULE) - 1, -1, -1):
        width = SCHEDULE[i]
        reserve = len(SCHEDULE) - i - 1
        setup = 44 if width == 24 else 37
        step = 61 if width == 24 else 23
        if reserve == 0:
            step += cube_ops + 43
        elif (i + 1) % len(WIDTHS) == 0:
            step += cube_ops + 24
        accepted = [value + step for value in continuation]
        worst = [0] * (LIMIT + 1)
        for (counts, _), effective in zip(rows[width], caps[i], strict=True):
            usage = [0] * (LIMIT + 1)
            for cost, count in counts:
                for cap in range(cost + reserve, LIMIT + 1):
                    usage[cap] += count * accepted[cap - cost]
            for remaining, cap in enumerate(effective):
                worst[remaining] = max(worst[remaining], usage[cap])
        continuation = [setup + value for value in worst]
    return continuation


def accelerated_work(caps):
    by_cap = state_work(caps)
    cube_ops = sum(operation_counts().values())
    assert len(BUDGETS) == RADIUS + 1
    return 1000 + COMPONENTS * sum(
        comb(M, d) * (38 + 54 + cube_ops + by_cap[cap])
        for d, cap in enumerate(BUDGETS)
    )


def counter_samples():
    from opmodel import Machine
    from construction import table_word

    assert tuple(OFFSETS) == (0, 13, 26, 39, 50)
    assert MEMORY_END < PAIR_END <= 1 << 32
    source, dest = ARRAYS["X0"], ARRAYS["C0"]
    program = tuple(cube_schedule(source, dest))
    charged = operation_counts()
    assert len({name for op, name, *_ in program if op != "store"}) + 4 <= 64
    rng = random.Random(910372)
    sample = (0, 1, MASK, *(1 << i for i in OFFSETS[:-1]),
              *(rng.getrandbits(M) for _ in range(8)))
    for value in sample:
        machine = Machine([0] * 2000)
        machine.mem[source] = value
        regs = {}
        for op, dst, *args in program:
            operands = [regs[a] if isinstance(a, str) else a for a in args]
            if op == "move":
                regs[dst] = operands[0]
            elif op == "store":
                machine.store(*operands)
            else:
                if op == "load" and operands[0] >= MEMORY_END:
                    machine.mem[operands[0]] = table_word(operands[0])
                regs[dst] = getattr(machine, op)(*operands)
        assert machine.mem[dest] == cube(value)
        assert machine.ops == sum(charged.values())
    return {"cube_ops": sum(charged.values()), "cube_op_counts": charged,
            "cube_replays": len(sample), "pair_memory_end": PAIR_END}
