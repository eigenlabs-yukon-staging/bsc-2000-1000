"""Twelve fixed chunks with one exact three-word cube data table per pair."""
from bisect import bisect_right
from collections import Counter
from itertools import accumulate, combinations
from math import comb
import random

from construction import ANCHOR_BITS, ARRAYS, BUDGETS, LIMIT, M, MASK, MASK64, MEMORY_END, RADIUS, WIDTHS, WORDS, cube


CHUNKS = (12,) * 11 + (11,)
OFFSETS = tuple(accumulate(CHUNKS, initial=0))
PAIRS = tuple(combinations(range(len(CHUNKS)), 2))
UNARY_OWNER = (1, 0) + (0,) * (len(CHUNKS) - 2)
BASES = tuple(accumulate((WORDS << (CHUNKS[i] + CHUNKS[j]) for i, j in PAIRS),
                         initial=MEMORY_END))
PAIR_END = BASES[-1]


def pair_word(address):
    if not MEMORY_END <= address < PAIR_END:
        raise ValueError("outside pair tables")
    index = bisect_right(BASES, address) - 1
    i, j = PAIRS[index]
    key, word = divmod(address - BASES[index], WORDS)
    a = (key & ((1 << CHUNKS[i]) - 1)) << OFFSETS[i]
    b = (key >> CHUNKS[i]) << OFFSETS[j]
    value = cube(a ^ b)
    if UNARY_OWNER[i] != j:
        value ^= cube(a)
    if UNARY_OWNER[j] != i:
        value ^= cube(b)
    return (value >> (64 * word)) & MASK64


def cube_schedule(source, dest):
    src = tuple(f"x{j}" for j in range(WORDS))
    for j, reg in enumerate(src):
        yield "load", reg, source + j
    lows = tuple(f"low{i}" for i in range(len(CHUNKS)-1))
    highs = {(j,CHUNKS[i]): f"high{j}_{CHUNKS[i]}" for i,j in PAIRS}
    for i, (start, width) in enumerate(zip(OFFSETS[:-1], CHUNKS, strict=True)):
        word, offset = divmod(start,64)
        if offset:
            yield "shr", "part", src[word], offset
        else:
            yield "move", "part", src[word]
        if offset+width>64:
            yield "shl", "t", src[word+1], 64-offset
            yield "or_", "part", "part", "t"
        if offset+width!=64 and start+width!=M:
            yield "and_", "part", "part", (1<<width)-1
        yield "mul", "scaled", "part", WORDS
        if i<len(lows):
            yield "move", lows[i], "scaled"
        for (j,k),reg in highs.items():
            if j==i:
                yield "shl", reg, "scaled", k
    acc=tuple(f"a{j}" for j in range(WORDS))
    for index, (i,j) in enumerate(PAIRS):
        yield "add", "key", lows[i], highs[j,CHUNKS[i]]
        yield "add", "address", "key", BASES[index]
        for word,reg in enumerate(acc):
            if word:
                yield "add", "address", "address", 1
            yield "load", "value", "address"
            if index:
                yield "xor", reg, reg, "value"
            else:
                yield "move", reg, "value"
    for j,reg in enumerate(acc):
        yield "store", None, dest+j, reg


def operation_counts():
    return Counter(op for op,*_ in cube_schedule(ARRAYS["X0"],ARRAYS["C0"]) if op!="move")


def state_work(caps):
    from dynamics import SCHEDULE, statistics
    bound = [0]*(LIMIT+1)
    costs = statistics(10)
    ops = sum(operation_counts().values())
    for i in range(len(SCHEDULE)-1,-1,-1):
        width = SCHEDULE[i]
        reserve = len(SCHEDULE)-1-i
        setup = 53 if width==24 else 46
        step = 61 if width==24 else 23
        if reserve==0:
            step += ops+121
        elif (i+1)%len(WIDTHS)==0:
            step += ops+60
        continuation = [w+step for w in bound]
        worst = [0]*(LIMIT+1)
        for (counts,_),effective in zip(costs[width],caps[i],strict=True):
            usage = [0]*(LIMIT+1)
            for price,count in counts:
                for cap in range(price+reserve,LIMIT+1):
                    usage[cap] += count*continuation[cap-price]
            for r,cap in enumerate(effective):
                worst[r] = max(worst[r],usage[cap])
        bound = [setup+w for w in worst]
    return bound


def accelerated_work(caps):
    bound=state_work(caps)
    cube_ops=sum(operation_counts().values())
    assert len(BUDGETS)==RADIUS+1
    return 1000 + sum(comb(ANCHOR_BITS,d)*(44+144+cube_ops+bound[r])
                      for d,r in enumerate(BUDGETS))


def counter_samples():
    from construction import table_word
    from opmodel import Machine

    assert WORDS==3 and M==143 and OFFSETS[-1]==M
    assert MEMORY_END<PAIR_END<=1<<32
    source,dest=ARRAYS["X0"],ARRAYS["C0"]
    program=tuple(cube_schedule(source,dest))
    cnt=operation_counts()
    assert len({name for op,name,*_ in program if op!="store"})+4<=64
    rng=random.Random(172125)
    samples=(0,1,MASK,*(1<<i for i in OFFSETS[:-1]),
             *(rng.getrandbits(M) for _ in range(8)))
    for sample in samples:
        mi=Machine([0]*2000)
        for j in range(WORDS):
            mi.mem[source+j]=(sample>>(64*j))&MASK64
        reg={}
        for op,dst,*args in program:
            arg=[reg[a] if isinstance(a,str) else a for a in args]
            if op=="move":
                reg[dst]=arg[0]
            elif op=="store":
                mi.store(*arg)
            else:
                if op=="load" and arg[0]>=MEMORY_END:
                    mi.mem[arg[0]]=table_word(arg[0])
                reg[dst]=getattr(mi,op)(*arg)
        result=sum(mi.mem[dest+j]<<(64*j) for j in range(WORDS))
        assert result==cube(sample)
        assert mi.ops==sum(cnt.values())
    return {"cube_ops":sum(cnt.values()),"cube_op_counts":cnt,
            "cube_replays":len(samples),"pair_memory_end":PAIR_END}
