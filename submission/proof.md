# Explicit [48,24]/[32,16] product code and its exact decoder

All bit indices are zero based, in least-significant-bit-first order. The
literal matrices, tables, decoder and both scores below concern the same
fully specified f and g. No stochastic ensemble, assumed weight distribution,
Monte Carlo estimate or external theorem supplies any premise.

## 1. Encoder

Let P48 and P32 be the literal lists in construction.py, of lengths 24 and
16 respectively. Each entry is a k-bit integer for its respective k. For
a k-bit block a and the corresponding list P, set

    parity(a) = XOR of P[i] for all set bits i of a,
    c(a) = a OR (parity(a) shifted left by k).

This is a completely determined linear injective map of k bits to 2k bits:
the first k transmitted bits equal a. The matrices are authoritative code
data. P48 originated as the low 24 bits of SHA-256('bsc48-v1-i'), read little
endian, for i=0,...,23. Its particular realized performance is enumerated,
not inferred from the seed or any cryptographic assumption. P32 originated
from systematic row reduction of binary monomial evaluation vectors of
degree at most two on five bits; that historical construction is not needed
for validity because the literal systematic matrix fully specifies the code.
No claimed minimum distance or external coding theorem is used.

Split the first 984 message bits into 41 blocks of 24, and take the remaining
16 as a final block. Apply c with P48 to the first 41 and P32 to the last.
Concatenate in that order. Length is 41*48+32=2000, message size is
41*24+16=1000, and blockwise systematic coordinates prove injectivity. Thus
f carries exactly 2^1000 messages. encode implements precisely this layout.

## 2. Every syndrome leader and the decoder

For either block type, write k=len(P). Define the 2k-column binary syndrome
map H by h[i]=P[i] for i<k, and h[k+j]=1<<j for 0<=j<k. Then

    H(x) = parity(x[0:k]) XOR x[k:2k],  and H(c(a))=0.

The last k columns are a basis, so the graph of all 2^k syndromes with edges
s -> s XOR h[i] is connected. leaders(P) in construction.py defines a unique
chosen leader for each syndrome by breadth-first search. Start at zero with
distance zero and correction zero; all other distances start at marker 255.
Process the FIFO queue in order, and neighbours in increasing bit order
0,...,2k-1. At first discovery of a neighbour, record distance[parent]+1,
and the parent's information correction XOR 1<<i if i<k, otherwise unchanged.
First discovery wins ties. Stop once all 2^k syndromes are discovered.

This terminates: each vertex is appended once, and the basis edges ensure
every vertex is reachable in at most k steps, below marker 255. The finite
queue cannot be exhausted before all reachable vertices have been found.
BFS processes vertices in nondecreasing distance; induction by layers proves
its assigned distances are shortest path lengths. A weight-w noise vector
supplies a w-edge path to its syndrome. Conversely XORing coordinate vectors
along a shortest path supplies a noise vector of that syndrome and weight
at most its length. A repeated coordinate could be cancelled and contradict
shortest length. Therefore assigned distances are exactly minimum coset
weights, with one deterministic chosen vector e_s for each syndrome s.

The recorded correction[s] is e_s[0:k], by induction on discovery parents.
The remaining half is uniquely determined by

    e_s[k:2k] = s XOR parity(correction[s]).

The full vector and tie rule are therefore determined even though only its
information half is stored. Define the block decoder d(x) to be

    (x XOR correction[H(x)]) AND (2^k-1).

Apply d separately to the 42 received blocks and concatenate their outputs
in the message layout above. This defines g on every 2000-bit input. It
uses neither randomness nor an iteration cap nor knowledge of sent data.

## 3. Exact error and claim 0

Fix any transmitted c(a), and write its received block as c(a) XOR z. Since
H(c(a))=0, its syndrome is H(z). The output equals a if and only if the
information half of z equals correction[H(z)]. For a given syndrome and
information half the parity half is unique, by the displayed equation.
Thus decoding succeeds if and only if z=e_[H(z)]. There is exactly one
successful noise pattern per syndrome. The chosen minimum-weight leader is
ML on BSC(1/10), but the calculation only needs this specific decoder.

Let A48[w] and A32[w] count chosen leaders of each minimum weight. The declared
experiment constructs ALL 2^24 and 2^16 leaders from the literal matrices and
counts their distances. Each histogram sums to its complete syndrome count.
Under BSC(1/10), any particular noise vector of weight w in length n has
probability 9^(n-w)/10^n. Consequently the exact per-block success values are

    s48 = sum_w A48[w]*9^(48-w) / 10^48,
    s32 = sum_w A32[w]*9^(32-w) / 10^32.

This counts precisely one successful pattern per syndrome, including tied
cosets only once. Every message has these same block success probabilities.
The 42 noise blocks are independent, so complete-message success is
s48^41*s32. Averaging over uniform U preserves it, and

    eps(f,g) = 1-s48^41*s32.

The submitted exact rational is the outward rounding

    S = floor(10^100*s48^41*s32),  B=1-S/10^100,

reduced in submission.json. Since S/10^100<=s48^41*s32, eps(f,g)<=B. verify.py
computes both entire histograms, the fractions s48,s32, S and B with integer
and Fraction arithmetic. It checks the computed B equals the manifest score
and that eps<=B<1. The denominator fits the verifier's rational digit limit.
Approximate decimal output is explanatory and is not evidence for the bound.

## 4. Complete fixed tables and memory

Input words 0..31 and output words 32..47 follow bsc-ops-v1 exactly. Fixed
read-only data begins at address 48. The layout, all in construction.py, is:

* T48=48, followed by 2^24 information-correction words for P48;
* T32=48+2^24, followed by 2^16 correction words for P32;
* LUT48=T32+2^16, followed by six consecutive 256-word byte-syndrome tables;
* LUT32=LUT48+6*256, followed by four such 256-word tables.

For a block's byte chunk j, table entry v is XOR of h[8*j+b] over set bits
b of the byte v. This determines every LUT entry exactly. XORing these
lookups yields the complete syndrome, by linearity of H. Each correction
word is the previously specified BFS information leader. All these tables
are built in full before decoding begins, as allowed by free table preload.

The highest table address plus one is END=16845360=48+2^24+2^16+10*256,
well below 2^32. There is no decoder working-memory allocation. Only registers
r0,r1,r2 are used. The host BFS arrays and queue construct the fixed tables;
they are not additional RAM state during decoding.

TableMemory represents already preloaded words with compact arrays and a
dictionary for input/output. Its get method returns a precomputed array
entry, without lazy search or arithmetic to generate its value. The counter's
Machine.load is unchanged and charges one load. This changes only the host
representation of RAM, not the RAM program or its primitive set. opmodel.py
is byte-identical to reference/opmodel.py; the experiment asserts the protected
SHA-256 9b440361c3fede634b8021705f6afd62e964c1d996045c91a16228acee296903.

## 5. Literal straight-line RAM program, correctness and claim 1

program() in decoder.py generates a literal instruction-tuple list before
reading any received data. A register operand is a one-tuple and a constant
operand is a literal integer. Each non-move tuple is exactly one listed
bsc-ops-v1 primitive. The only moves assign literal zero, free by the rules.
The list has no runtime branch, loop, call or return instruction. All Python
loops, index calculations, branches and calls in the generator are compile
time operations. The host interpreter decode executes that RAM instruction
list on the unchanged reference counter. Its host dispatch is not a second
RAM program; every operation on a received value is an emitted counted
primitive. This is an explicit finite program, not a cost estimate.

For each literal block position the generator emits:

1. Extract its 48 or 32 input bits using a load, optional right shift, and
   when crossed a second load, left shift and OR. AND with the n-bit mask.
   Cost is 2 plus one for a nonzero starting offset plus three if crossed.
2. Free move zero to r2. For each literal byte emit shift-right, AND 255,
   add literal LUT base, load, XOR into r2: five primitives per byte.
3. Add literal correction-table base to r2, load the correction, XOR with
   the received block and AND the k-bit payload: four primitives.
4. Place the k decoded bits at its literal output position by load,
   shift-left, OR and store. For a crossing add load, shift-right, OR,
   store for the next word: four plus four if crossed.

Masked extraction produces exactly the received block even across input
word boundaries. Output memory begins zero. Output block bit intervals are
disjoint; the first shift truncates to the lower word and the optional second
store supplies precisely the remaining high part. OR placement therefore
builds the exact concatenated message. No operation sets output bits above
999, so required padding remains zero. Hence this instruction list computes
exactly g on all received words, including arbitrary large-noise patterns.

There are 42 blocks and 41*6+4=250 byte chunks. The fixed literal schedule
has 31 nonzero input offsets, 21 input crossings and 10
output crossings, counted from positions (48*b,24*b), b=0,...,40, and the
last block's position (1968,984). Thus exact primitive count is

    42*(2+4+4) + 250*5 + 31 + 3*21 + 4*10 = 1804.

This same instruction list executes on every input. Consequently W=1804
is the worst-case count, not a sampled average. For the claimed
L=10817/1000, the positive integer comparison

    1804^1000 <= 2^10817

is exactly equivalent to W<=2^L. verify.py checks it using integers and counts
the entire emitted instruction list independently of received values. It
also runs that list on zero, one, seeded arbitrary received words and encoded
messages; all sample counts and outputs must agree. Those samples illustrate
the proof; the fixed-list argument proves the all-input worst case.

## 6. Evidence, provenance and limits

The declared verify run is standard-library Python only. It builds every
leader, counts the code-specific histograms, computes the rational score,
checks every byte LUT and runs the exact decoder on the byte-identical counter.
The claim depends on its captured completed run within the protected
120-second/1024-MiB sandbox limits. It does not rely on an unrun computation.

The larger [50,25] exploratory table exceeded the runtime limit and is not
used in either claim. An exploratory QR48 alternative was also checked; the
literal seeded matrix here gives the stronger tested success probability.
The prior first-copy entry is a separate construction; no result from it is
an accepted premise for this proof. Both claims are proved here from the
current explicit code. This is still a weak full-block code, far from the
roughly 0.0585 nonexplicit random-coding reference. Structural checking,
sandbox execution and advisory critique are distinct from hosted review
and human confirmation.
