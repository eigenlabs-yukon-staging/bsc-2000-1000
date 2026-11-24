# A 60-bit systematic product with bitset-certified decoding success

## Encoder and syndrome

Bits are numbered from zero, least significant first. Split the message
into 33 consecutive 30-bit integers and a final 10-bit integer. For a
width k component define the linear map P by its k width-k columns.
For k=30, column i is the low 30 bits of the little-endian integer
given by the first four bytes of SHA-256 of the ASCII string
`wide-product-a:30:` followed immediately by the unpadded decimal integer
i, for i=0 through 29. The hash fixes literal data; no statistical
property of the hash is assumed. For k=10 the columns in increasing
index order are (980,731,346,811,541,494,407,357,441,744).
P(a) XORs exactly those columns for which a has a one bit.
The encoder sends each a followed by P(a) and concatenates in component
order. There are 33*30+10=1000 information bits and 2000 transmitted
bits. The systematic bits uniquely recover the message without noise,
so there are exactly 2^1000 encoded words.

For a received component written as (a,v), its syndrome is
S(a,v)=P(a) XOR v. The kernel of S consists precisely of the encoded
components. Its 2k columns, in coordinate order, are the k columns of
P followed by the k unit vectors. For any syndrome s, choose the
least-weight vector e of syndrome s. Among equally weighted vectors,
use the lexicographically first increasing tuple of flipped coordinates.
Denote its systematic portion by C(s). A candidate always exists, since
a vector supported on the last k coordinates can realize any syndrome.
In particular the following offline procedure constructs every constant
table exactly. Fill the P table by direct column XOR with the fixed high-bit marker specified below. Initialize all C
entries as unassigned, enumerate increasing tuples by weight from 0
through k, and set C at the syndrome of each tuple on first encounter
to the low k bits of its bit vector. Ignore subsequent encounters.
Enumeration halts after a finite number of steps, and all entries are
assigned. This is a definition of preloaded data; verification of the
error bound does not enumerate the complete correction table.

For each received (a,v) the decoder outputs a XOR C(S(a,v)). This is
the systematic part after removal of the chosen vector. The corrected
vector has syndrome zero and is consequently an encoded component.
The full decoder concatenates these outputs.

## Error bound

Since the sent component has zero syndrome and S is linear, the syndrome
of a received component equals the syndrome of the actual noise. The
decoded component is correct exactly when the chosen and actual noise
vectors are equal. If the output is correct, the corrected encoded word
must be the transmitted word by systematicity, hence the XOR of the
received and corrected words is precisely the actual noise. Conversely,
removing the actual noise recovers the sent word.

Let R_t be the set of syndromes of vectors of weight at most t. Then
R_0={0}, and R_(t+1) = R_t union union_j (R_t XOR column_j).
Indeed, toggling one coordinate increases weight by at most one.
Every positive-weight vector of weight at most t+1 can also be formed
by toggling a nonzero coordinate of a vector of weight at most t.
By induction, |R_t|-|R_(t-1)| counts precisely the chosen leaders of weight
t, with the weight-zero count being one.

The declared error experiment changes basis for the
wide syndrome using the thirty fixed row masks in
BASIS_ROWS. A transformed bit is the parity of the
old syndrome intersected with its row mask. Exact
binary elimination verifies full rank. Thus the map
permutes all syndromes, and applying it to all sixty
columns preserves each reachability count. It is a
host counting optimization and has no effect on the
counted decoder.

To initialize through weight five, the verifier
enumerates all increasing tuples of coordinate indices
at each weight, XORing their transformed columns
and setting the corresponding bit. Through weight ten,
a step from exact R_(t-1) splits the sixty columns
into a selected list and its complementary list by
their literal order after sorting by low bit count.
The selected list has length 28,32,36,40,60 in the
five steps. Start with a copy of R_(t-1). Apply XOR
translations by selected columns to that preceding
set and union the results. This reaches every syndrome
of any pattern of weight t using a selected coordinate,
by removing one such coordinate. Then enumerate
the weight-t tuples using only complementary columns,
which reaches all remaining patterns. Every included
syndrome is realized by an error of weight at most t.
Therefore these steps compute exactly R_t.

The low thirteen transformed syndrome bits are
stored as bit positions of each integer bitset row;
the high seventeen bits are its row index. XOR with
a column permutes the rows by its high bits and
permutes bit positions by its low bits. For low
index bit j use step=2^j. The mask
((2^(2^13))-1)/(2^step+1) selects alternating
lanes of width step, since 2^13 is a multiple of
2*step. Its shift left by step selects complementary
lanes. Shift the masked values in opposite directions
by step and OR to exchange lanes, which toggles
that index bit. Composing these exchanges realizes
the column's low XOR permutation. The bitset union
reads only the preceding rows. For an at-least
three-quarters full preceding set, it runs by
destination row so union can stop for that row
once all bits are set. XOR is an involution, so
this gives the same result as source-row traversal.
Complementary enumeration sets bits in the new
rows after the translations.

The tail counts exact reachability by enumerating
all tuples at successive weights to saturation.
Both histograms come from differences of cumulative
counts. The wide counts through weight ten are
(1,60,1770,34220,487565,5447068,48692266,308302482,
645060687,65715621,84), summing to 2^30. The tail
counts are (1,20,190,740,71,2), summing to 2^10.
The verifier checks both saturations.

For a weight-w error in length n the channel
probability is 9^(n-w)/10^n. Let A be the weighted
sum of wide leader counts with weights 9^(60-w),
and T that of the tail counts with weights
9^(20-w). Each component's exact success is its
weighted sum divided by the relevant power of ten.
The channel flips in separate components are
independent for every message. Thus whole-message
success equals A^33*T/10^2000. The exact arithmetic
checks

    1-A^33*T/10^2000 <= 499976397292437707/500000000000000000.

## RAM implementation and cost

For width k the fixed memory contains an address table
with each entry containing C_base OR P(a) and a correction
table with each entry containing C(s). For k=30 the two
bases are P_base=2^31 and C_base=2^30. For k=10 the
bases are P_base=2^11 and C_base=2^10. Each range has
2^k words. They are disjoint, all start at addresses at
least 48, and the largest exclusive end is 3*2^30 below
2^32. Each word fits the unsigned register width. The
fixed C_base marker has no bits in the k low positions
because C_base=2^k in both cases.

The RAM program is the finite tuple sequence emitted by schedule() in
construction.py. Every tuple other than move denotes one named RAM
primitive, and each variable name denotes a register. Operands given
as integers are literal constants. An ordinary tuple assigns the
primitive result to the destination register. Store has no destination
register. Move copies a register or literal for free. All loops and
conditions in schedule() run at code generation time, independent of
any received values, and produce the same fixed sequence. The runtime
has no branches, calls or returns. It halts after the stores.

There are 32 input registers Y, 16 output accumulators U, and four
temporaries R, T, A, P, for 52 registers. The initial loads and free
zero moves establish the input registers and output accumulators.
For component starting at an offset of s within word Y[j], R first
holds Y[j] shifted right by s (or copied for s=0), with Y[j+1]
shifted left by 64-s and ORed into R if the component crosses a word.
The two terms in this OR have disjoint positions. The low 2k bits of
R are exactly the component's received bits. The program masks the
low k bits into A and shifts R right by k to start the received parity
in P. It masks P to k bits unless R already has exactly 2k valid bit
positions after shifting and required zero filling. In either case the
value in P now equals the received v. One address addition
and load obtain the entry C_base OR P(A). XOR this with
v. Since both v and P(A) have k bits and C_base=2^k,
the result is exactly C_base+S. Load directly at that
address to obtain C(S), and XOR with A to decode.
Both loads and the single address addition are counted.
Every computed address stays within its own table range.

The components are processed in order. For a component with output
offset s inside a word, shifting the k-bit result left by s and ORing
it into the accumulator writes its low portion. If s is zero, it can
instead replace the untouched zero accumulator. If a component crosses
the boundary, right shifting by 64-s initializes the next accumulator
with its high portion. The ascending disjoint output intervals imply
that both direct assignments initialize words with no earlier output.
Final stores write the 16 output words and leave the last 24 padding
bits zero. Hence this program computes exactly the analyzed decoder
for all inputs.

The emitted sequence has 100 loads, 66 ANDs, 79 right shifts, 34
additions, 68 XORs, 61 left shifts, 61 ORs and 16 stores. The exact
total is 485 on every input. The declared ops experiment counts that
list and verifies 485^500 <= 2^4461. For selected sample inputs, with at
most two flips per component, it precomputes the exact values at all
accessed table addresses. Exhaustive enumeration through weight two
suffices to identify the chosen leader of each such syndrome. Extra
addresses are unused by those executions and the test raises an error
if they are accessed. The samples exercise the same primitive sequence
through the byte-identical reference counter. The sample preload is
only a test fixture; the mathematical definition above fixes the full
preloaded tables and the input-independent count does not rest on
sampling.
