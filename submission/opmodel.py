"""bsc-ops-v1: the organizer's reference counter for the decoding-ops objective.

The cost model is defined in rules.md ("Decoding-ops model"); this file is an
executable restatement of it. A submission may copy this file byte-for-byte
into submission/ and run its decoder on it in a declared experiment. One run
counts the operations the decoder executes on one received word. It is
evidence for the worst-case bound, never the bound itself, unless the proof
shows the decoder's control flow (and so its op count) is the same on every
received word. Python 3 standard library only.

Memory layout (64-bit words, addresses 0 .. 2^32 - 1, initially zero):
  0..31    received word y: bit i of y is bit (i % 64) of word i // 64 (LSB first);
           bits 2000..2047 are zero. Preloaded at no cost.
  32..47   output u: bit i of u is bit (i % 64) of word 32 + i // 64; bits
           1000..1023 must be zero. Written with store, which is counted.
  48..     the decoder's fixed read-only tables (preloaded at no cost, part of
           the decoder's specification) and its working memory.
"""

MASK = (1 << 64) - 1
MEMORY_WORDS = 1 << 32
INPUT_BASE, OUTPUT_BASE, TABLE_BASE = 0, 32, 48
N, K = 2000, 1000


class OpLimit(Exception):
    pass


class Machine:
    """Every method below is one primitive and costs exactly 1 op.

    Values are unsigned 64-bit words (Python ints in 0..2^64-1) or binary64
    floats. Constants written in the program and moves between the 64 registers
    are free; the caller keeps register values in local variables and must
    route every data-indexed access through load/store.
    """

    def __init__(self, received_bits, tables=None, rand_bits=None, op_limit=None):
        if len(received_bits) != N or any(b not in (0, 1) for b in received_bits):
            raise ValueError("received word must be 2000 bits")
        self.ops = 0
        self.op_limit = op_limit
        self.mem = {}
        for i, b in enumerate(received_bits):
            if b:
                self.mem[INPUT_BASE + i // 64] = self.mem.get(INPUT_BASE + i // 64, 0) | (1 << (i % 64))
        for addr, word in (tables or {}).items():
            if not TABLE_BASE <= addr < MEMORY_WORDS:
                raise ValueError("tables live at addresses 48 .. 2^32 - 1")
            self.mem[addr] = word & MASK
        self._rand = rand_bits

    def _tick(self):
        self.ops += 1
        if self.op_limit is not None and self.ops > self.op_limit:
            raise OpLimit(self.ops)

    @staticmethod
    def _addr(a):
        if not 0 <= a < MEMORY_WORDS:
            raise IndexError("address outside 0 .. 2^32 - 1")
        return a

    # memory
    def load(self, a):
        self._tick(); return self.mem.get(self._addr(a), 0)

    def store(self, a, v):
        self._tick(); self.mem[self._addr(a)] = v & MASK if isinstance(v, int) else v

    # integer word operations (results reduced mod 2^64)
    def add(self, a, b): self._tick(); return (a + b) & MASK
    def sub(self, a, b): self._tick(); return (a - b) & MASK
    def mul(self, a, b): self._tick(); return (a * b) & MASK
    def divu(self, a, b): self._tick(); return a // b
    def modu(self, a, b): self._tick(); return a % b
    def and_(self, a, b): self._tick(); return a & b
    def or_(self, a, b): self._tick(); return a | b
    def xor(self, a, b): self._tick(); return a ^ b
    def not_(self, a): self._tick(); return ~a & MASK
    def shl(self, a, s): self._tick(); return (a << s) & MASK if s < 64 else 0
    def shr(self, a, s): self._tick(); return a >> s if s < 64 else 0
    def popcount(self, a): self._tick(); return bin(a).count("1")
    def ctz(self, a): self._tick(); return (a & -a).bit_length() - 1 if a else 64
    def eq(self, a, b): self._tick(); return int(a == b)
    def ltu(self, a, b): self._tick(); return int(a < b)
    def lts(self, a, b):
        self._tick()
        sa = a - (1 << 64) if a >> 63 else a
        sb = b - (1 << 64) if b >> 63 else b
        return int(sa < sb)
    def select(self, c, a, b): self._tick(); return a if c else b

    # binary64 operations (IEEE round-to-nearest, as Python floats)
    def fadd(self, a, b): self._tick(); return float(a) + float(b)
    def fsub(self, a, b): self._tick(); return float(a) - float(b)
    def fmul(self, a, b): self._tick(); return float(a) * float(b)
    def fdiv(self, a, b): self._tick(); return float(a) / float(b)
    def flt(self, a, b): self._tick(); return int(float(a) < float(b))
    def fmin(self, a, b): self._tick(); return min(float(a), float(b))
    def fmax(self, a, b): self._tick(); return max(float(a), float(b))
    def fabs(self, a): self._tick(); return abs(float(a))
    def itof(self, a): self._tick(); return float(a)
    def ftoi(self, a): self._tick(); return int(a) & MASK

    # control flow: call once for every branch or jump taken or not taken,
    # including each loop back-edge and each exit test, and for call/return.
    def branch(self, taken=True): self._tick(); return taken

    # randomness: 64 fresh uniform bits
    def rand(self):
        self._tick()
        return next(self._rand) & MASK if self._rand is not None else 0

    def output_bits(self):
        """The decoded message, read for free after the decoder halts."""
        if self.mem.get(OUTPUT_BASE + 15, 0) >> 40:
            raise ValueError("output bits 1000..1023 must be zero")
        return [(self.mem.get(OUTPUT_BASE + i // 64, 0) >> (i % 64)) & 1 for i in range(K)]
