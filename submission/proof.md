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
table exactly. Fill the P table by direct column XOR. Initialize all C
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

The declared error experiment implements this recurrence using integers
as exact bitsets. For low group width l, a row is indexed by the high
k-l syndrome bits and holds 2^l low-index positions. XOR with a column
permutes rows by its high bits. XOR with bit j of the column's low
part exchanges each adjacent pair of lanes of width step=2^j inside a
row. The mask ((2^(2^l))-1)/(2^step+1) is the sum of lane masks
(2^step-1) shifted by multiples of 2*step. The division is exact
because 2^l is a multiple of 2*step. Its left shift by step selects
the complementary lanes. Masking and shifting the two sets of lanes
in opposite directions and ORing them implements that exchange.
Composing for each set bit in the low column XORs the entire low
index. OR into a fresh row array initially copied from R_t computes
the union, while reads use the unchanged preceding row array.
When at least three quarters of the states are already present, the
experiment instead traverses by destination row, using the same XOR
translation in the reverse direction, which equals the forward
translation because XOR is its own inverse. In this traversal a row
can stop receiving updates once all its bits are set.
Population counts followed by differences give the histogram. The low
group widths are 13 and 5 respectively for the 30- and 10-bit
syndromes. All arithmetic here is integer arithmetic.

Writing h_w for the wide histogram and t_w for the tail histogram,
the experiment computes h_0 through h_9 and all t_w through
saturation. The wide counts are
(1,60,1770,34220,487565,5447068,48692266,308302482,645060687,65715621).
They sum to 1073741740, leaving 84 syndromes out of this success
certificate. The tail counts are (1,20,190,740,71,2) and sum to 1024.
For a particular weight-w error vector in length n the channel
probability is 9^(n-w)/10^n. Therefore with
A = sum_(w=0..9) h_w * 9^(60-w) and
T = sum_w t_w * 9^(20-w), per-component success probabilities are at
least A/10^60 and exactly T/10^20 respectively. The wide inequality
uses only the certified leaders; the defined decoder is complete on
every syndrome. For every fixed message, the component success events
depend on disjoint independent sets of channel flips. Thus whole
message success is at least A^33*T/10^2000 for every message and
also on average. The exact verifier arithmetic checks
1-A^33*T/10^2000 <= 499976397292483003/500000000000000000.

## RAM implementation and cost

The fixed memory contains two disjoint direct tables for each width,
one containing P(a) for every a and one containing C(s) for every s.
The bases appear in table_addresses() in construction.py, starting at
48, in the order P30, C30, P10, C10, with range lengths 2^30, 2^30,
2^10, 2^10. The last address plus one is 2147485744, below 2^32.
Each table value occupies one word and fits in k bits.

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
value in P now equals the received v. A table address addition and
load compute P(A); XOR with v computes S. A second address addition
and load compute C(S); XOR with A yields the decoded component.
All dynamic addresses fall in the fixed ranges and both additions and
both table loads are counted.

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

The emitted sequence has 100 loads, 66 ANDs, 79 right shifts, 68
additions, 68 XORs, 61 left shifts, 61 ORs and 16 stores. The exact
total is 519 on every input. The declared ops experiment counts that
list and verifies 519^50 <= 2^451. For selected sample inputs, with at
most two flips per component, it precomputes the exact values at all
accessed table addresses. Exhaustive enumeration through weight two
suffices to identify the chosen leader of each such syndrome. Extra
addresses are unused by those executions and the test raises an error
if they are accessed. The samples exercise the same primitive sequence
through the byte-identical reference counter. The sample preload is
only a test fixture; the mathematical definition above fixes the full
preloaded tables and the input-independent count does not rest on
sampling.
