# Exact first-copy point

## Complete code and decoder

Index bits from zero. Let f(u)[i]=u[i mod 1000] for 0<=i<2000.
Let g(y)[i]=y[i] for 0<=i<1000. Both functions are deterministic and
fully specified. The first 1000 coordinates make f injective, hence there
are exactly 2^1000 codewords, each of length 2000. The decoder returns
exactly 1000 bits and has no access to either the message or the noise.

## Claim 0: exact block error

For every fixed u, the event g(f(u) XOR Z)=u is precisely the event that
all of Z[0],...,Z[999] are zero. The last 1000 noise coordinates have
no effect on the decoder. Independent BSC(1/10) noise therefore gives
success probability (9/10)^1000 for every u, and averaging over uniform U
preserves this value. Thus eps(f,g)=(10^1000-9^1000)/10^1000 exactly.
To stay within the verifier's 256-digit rational format, claim 0 uses
an exact rational outward rounding rather than that 1001-digit fraction.
Set S=floor(9^1000*10^100/10^1000) and B=1-S/10^100, reduced as printed
in submission.json. Since S/10^100 <= (9/10)^1000, the exact error is at
most B. exact_check.py performs this integer division, reduces B with
Fraction, prints its numerator and denominator, and checks both
1-(9/10)^1000<=B and B<(10^46-1)/10^46. These integers match claim 0's
score. No floating point, ensemble assertion or empirical error estimate
enters the claim. The quotient definition determines B without ambiguity.

## Claim 1: complete RAM program and worst-case count

The machine and input/output layout are exactly bsc-ops-v1 in rules.md.
The following is a macro describing a fully unrolled straight-line program:

```
For each literal j=0,1,...,14, emit these two instructions:
    r := load(j)
    store(32+j, r)
r := load(15)
r := and(r, 1099511627775)  # literal 2^40-1
store(47, r)
HALT
```

The phrase 'for each literal' is code generation, not a runtime loop.
All 33 instructions have fixed addresses and literal constants. No
runtime indexing, branch, loop, call or return occurs. Moves, constants
and halt are free by the rules. There are 16 loads, 16 stores and one
64-bit and, so exactly W=33 operations on every received word. Only one
64-bit data register is necessary; even allowing a second address
register is well below 64. Input addresses 0..31 and output addresses
32..47 are the only addresses; there are no tables or working-memory
allocations. The 2^32-word memory limit is respected.

Input words 0..14 hold y[0:960]. Input word 15 holds y[960:1024]. The
mask keeps its low 40 bits, namely y[960:1000], and zeros the remaining
24 output padding bits. Thus this exact program computes g on every
received word, not merely on codewords or low-weight errors.

For L=12611/2500, W<=2^L is equivalent (both sides positive) to
33^2500<=2^12611. ops_check.py verifies this exact integer comparison
and checks 12611*20<101*2500. The latter proves our L is strictly below
the starter's 101/20. The counter file opmodel.py is copied byte-for-byte
from reference/opmodel.py with the published SHA-256. Its calls in
ops_check.py denote primitive instructions of the RAM program; the host
Python loop is the macro unroller, not additional RAM control flow, as
in the supplied starter. Test cases illustrate instruction counts and
correct outputs; the input-independent straight-line proof establishes
the worst case. Neither error nor cost claims rely on simulation.

## Provenance and scope

This entry deliberately keeps the official starter construction and its
33-instruction implementation. The contribution is a rigorously outward-rounded error score
and a tighter rational complexity certificate, not a new coding algorithm.
Both exact rational scores strictly improve the official starter point. Local
structural checking and advisory review are distinct from the configured
review and human confirmation required for official acceptance.
