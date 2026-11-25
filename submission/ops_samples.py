"""Sample the fixed 283-operation decoder with the protected RAM counter."""
from itertools import combinations
from collections import Counter
import random

from construction import (
    MAIN_COLUMNS, SYNDROME_COLUMNS, OFFSET, TAIL_BASE,
    MASK32, parity, encode, schedule, table_word, tail_leaders,
)
from opmodel import Machine

MASK64=(1<<64)-1


def short_main_leaders():
    out={}
    for weight in range(3):
        for locations in combinations(range(64),weight):
            s=0
            for j in locations:
                s^=SYNDROME_COLUMNS[j]
            if s not in out:
                out[s]=sum(1<<j for j in locations)
    return out


def bits(x,n):
    return [(x>>i)&1 for i in range(n)]


def analyze_word(y,program,main,tail):
    tables={}
    in_words=[(y>>(64*j))&MASK64 for j in range(32)]
    out=[0]*16

    def preload(key):
        value=table_word(key,main,tail)
        tables[key]=value
        return value

    def read(address):
        if address<32:
            return in_words[address]
        if address<48:
            return 0
        return preload(address)

    for j,r in enumerate(in_words[:31]):
        a=r&MASK32
        check=read(a)
        index=(check^r)>>32
        correction=read(index)
        output_word,lane=divmod(j,2)
        out[output_word]|=((a^correction)&MASK32)<<(32*lane)
    out[15]|=preload(TAIL_BASE+in_words[31])
    expected=sum(v<<(64*j) for j,v in enumerate(out))
    machine=Machine(bits(y,2000),tables=tables)
    reg={}
    for op,dst,*args in program:
        used=[reg[x] if isinstance(x,str) else x for x in args]
        if op=="move":
            reg[dst]=used[0]
        elif op=="store":
            machine.store(*used)
        else:
            if op=="load" and used[0]>=48:
                assert used[0] in tables
            reg[dst]=getattr(machine,op)(*used)
    assert machine.output_bits()==bits(expected,1000)
    assert machine.mem[47]>>40==0
    assert machine.ops==283


def check_runtime(exponent):
    program=tuple(schedule())
    counts=Counter(op for op,*_ in program if op!="move")
    assert sum(counts.values())==283
    assert 283**exponent.denominator<=2**exponent.numerator
    rng=random.Random(550987)
    examples=[encode(0),encode((1<<1000)-1)]
    for _ in range(8):
        now=encode(rng.getrandbits(1000))
        for j in range(31):
            for i in rng.sample(range(64),rng.randrange(3)):
                now^=1<<(64*j+i)
        now^=rng.getrandbits(16)<<1984
        examples.append(now)
    addresses=(0,3,31,32,47,48,TAIL_BASE-1,TAIL_BASE,TAIL_BASE+65535)
    outside=0
    safe_data=0
    for j in range(31):
        a=rng.randrange(48,TAIL_BASE)
        v=parity(a,MAIN_COLUMNS)^OFFSET^addresses[j%len(addresses)]
        outside|=(a|(v<<32))<<(64*j)
        b=addresses[j%len(addresses)]
        safe_data|=(b|(rng.getrandbits(32)<<32))<<(64*j)
    examples.append(outside|(65535<<1984))
    examples.append(safe_data)
    class PartialLeaders(dict):
        def __missing__(self, index):
            return ((index*0x44281317)^0xaafd27) & MASK32
    main,tail=PartialLeaders(short_main_leaders()),tail_leaders()
    for y in examples:
        analyze_word(y,program,main,tail)
    return {"ops":283,"counts":counts,"replayed_words":len(examples)}
