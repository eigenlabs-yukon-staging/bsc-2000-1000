# Word-aligned 64-bit components in 283 operations

## Code and local syndrome leaders

Use zero-based binary indices, low bit first. Divide
the first 992 message bits into 31 consecutive blocks
of 32 bits each. For each block a, send a followed
by its 32-bit parity P(a). To define P, start
with the two 16-bit rows in construction.py. The first character
is bit zero. Take their cyclic left rotations by i=0..15, packing
the rotation of the first row below the rotation of the second.
For columns 16..31 take the rotations of the transposed second
and transposed first rows in that order. Transposing here replaces
bit i of a first row by its bit (-i mod 16). This makes the four
circulant parity blocks explicit. The base row strings are from
code C64,1 in Anev, Harada and Yankov's Table 1; the score
proof recalculates the syndrome counts of the actual modified
code. [New extremal singly even self-dual codes of lengths 64 and 66](https://arxiv.org/abs/1708.05950)

Now XOR bit 27 of parity column i exactly when the ith bit
of the fixed data mask 0x2e27aa03 is one. This changes one parity
coordinate of each codeword by the corresponding linear form.
For any block a, XOR the resulting columns for its nonzero
data bits. This fixes a linear systematic [64,32] component.

Send the last eight message bits followed by the
eight-bit parity from the column list TAIL_COLUMNS.
The layout has length 31*64+16=2000, with exactly
1000 data coordinates and 2^1000 distinct messages.

A main received block (A,V) has syndrome
s=P(A) XOR V. The kernel of this linear map is
the code, and its 64 columns are the 32 P columns
and 32 unit bits in order. For each syndrome choose
a least-weight error, breaking ties by the first
lexicographic increasing tuple of flipped coordinates.
This is complete and explicit: enumerate those
tuples by weight through 32, assigning each newly
seen syndrome its first error. All syndromes occur
because the parity unit columns span them. Let
C(s) be the selected error's 32-bit data part.
Ordinary main decoding returns A XOR C(s).
The tail code likewise chooses minimum weight
errors by breadth-first search of its syndrome
graph, with edge order as in tail_leaders. The
shortest path length equals the minimum error
weight because toggling an index twice cancels
and every error tuple is a walk. These choices
fix all table words.

## Shared parity and correction table

Let rho=OFFSET, H=2^32-2^16, and C(s) be the selected
correction's 32 low data bits. Use one table at addresses
48 through H-1, inclusive. At address r its word contains
C(r XOR rho) in the low half and P(r) XOR rho in the
high half. Both fields are fully determined for each r.
They serve distinct lookups of the same bank. Addresses
H through 2^32-1 are unchanged tail-decoder entries
indexed by the sixteen received tail bits. Their values
contain only a decoded eight-bit field at bits 32..39.
The two fixed ranges are disjoint and lie within RAM.

For each received main word Y with low half A and high
half V, load Y, mask it to form A and read the word at
address A. XOR that word into Y and shift right by 32.
When A is in 48..H-1, the result is r=P(A) XOR V XOR rho.
Whatever the result, it is a valid 32-bit address.
Read its word and XOR it into the original Y. When
both addresses are in the shared table, the low half
is then A XOR C(P(A) XOR V), the desired ordinary
systematic correction.

An even component masks to the low half so the upper
32 bits of its output lane are zero. An odd component
shifts left by 32, retaining only its low half in the
other lane, and ORs the two lanes together. This holds
for all input words, including addresses outside the
shared bank. The input prefix reads received words
and the output prefix reads zeros until the final stores.
The tail uses four operations to read its received word,
form a table address, load the shifted decision and OR
it into the last output word. The tail read has a valid
sixteen-bit index because input padding is zero. The
result's upper 24 bits are zero for all inputs.

There are 31 main components of eight operations each,
15 additional ORs for odd components, four tail
operations and sixteen final stores. Thus W=283.
There are at most sixteen output registers and eight
other named registers. Branches and loops in schedule()
take place at generation time from known literal indices.
The emitted sequence is fixed. No loops, address work or
conditional jumps remain implicit at decoder runtime.

## Exact syndrome certificate

A linear codeword has syndrome zero, so the
received syndrome is exactly that of the noise.
Ordinary correction succeeds exactly when the
noise is the selected leader. Excluding main
correction keys below 48 or at least H therefore
removes exactly one success contribution per
excluded syndrome. For a syndrome of minimum
weight w whose key is safe, that contribution is
9^(64-w)/10^64. This probability holds for each
transmitted message, with the same chosen error.

The exact histogram computation uses the 32 row
masks BASIS in verify.py as a new binary basis.
Its ith output bit is the parity of the input
syndrome AND row i. The declared rank check
eliminates the fixed masks and verifies rank 32.
Thus this is a bijection on syndromes and its
linear action on columns preserves the weight
of a minimum error for each translated syndrome.
The top two transformed bits define four disjoint
sectors. 46 transformed columns have top
two bits zero and thus preserve all four sectors
under XOR. The other 18 are available
for direct enumeration but are never used in
a sector-local bitset translation.

For a sector let R_t be its transformed syndromes
attained by errors of weight at most t. Enumerate
all masks through weight five, keeping only
syndromes in that sector. For each weight t>=6
choose a subset E of the preserving columns
using the smallest low-bit counts in the exact
sort order of the checker. The counts selected are the smaller of the preserving
column count and 33, 36, 40, 44 at weights 6, 7, 8, 9.
At later weights all the preserving columns are selected. Begin with R_(t-1),
then OR its translations by all selected columns.
Every single selected flip preserves the sector.
If a weight-t error reached this sector and
uses a selected coordinate, removing that
coordinate gives an error of weight t-1 in the
same sector. Thus these translations reach
its syndrome. For weight-t errors with no
selected coordinate, explicitly enumerate all
t-tuples from the complementary columns and
keep those syndromes of the sector. Conversely
every reached syndrome has an error of weight
at most t by construction. Induction proves
the sets computed are exactly R_t.

Each sector stores 131072 exact integer bitset
rows, each with 2^13 bit positions. Low transformed bits index the bit
within the row, with the other seventeen
local bits indexing the row. XOR by a column's
high local part permutes rows. To XOR low
index bit j, use step=2^j and the alternating
lane mask ((2^(2^13))-1)/(2^step+1). This is
an integer since the row width is a multiple
of 2*step. One mask and its shift left by step
isolate complementary lanes. Shift their
masked values in opposite directions and OR
to exchange them. Composition realizes each
column's low permutation. The updates read
only the preceding exact rows. With enough
bits set they can traverse by destination
row and stop once a row is full, which is
equivalent because translation by XOR is
its own inverse. This is host exact
counting, separate from the RAM decoder.

Construct an exclusion mask in the same
representation by taking every actually
excluded address r, mapping r XOR rho
through the basis, and including it in its
transformed sector. Its intersections
with R_t count the excluded syndromes
already reached. Differences of consecutive
total and exclusion counts give respectively
all minimum weights and excluded minimum
weights. Both histograms are stored in
verify.py and independently asserted by
each run for its designated sector. The
four declared runs together cover the
entire syndrome space and all 48+2^16
excluded keys, without assuming an
unverified histogram. Each run verifies saturation at or before weight twelve.

Let A be the sum of 9^(64-w) over the
safe syndrome leaders counted in these
four histograms. Let T be the sum of
9^(16-w) over the tail leaders, which the
checker enumerates completely by BFS.
For fixed noise on one main component, its received
data half is the XOR of its uniform 32-bit message
and fixed systematic noise. It is therefore uniform
over all 2^32 values. This holds conditional on every
noise pattern, including any chosen syndrome leader
with safe correction address. Thus the probability of
a safe parity address conditional on such noise is
D/2^32, where D=2^32-48-2^16. The probability of
successful main decoding with both addresses safe is
at least A*D/(10^64*2^32). This lower bound counts
only chosen leaders with safe keys and safe data halves.

The components use disjoint input bits under the uniform
message and disjoint independent flips. On the intersection
of the certified local success events the packed output
equals the full message. The tail success probability
is T/10^16. Consequently block error is at most

    1-(A*D)^31*T/(2^992*10^2000).

The four declared computations check this inequality
against the exact rational claim. Sector zero counts
all 283 operations, checks W^L_den <= 2^L_num, and
runs the byte-identical reference model on selected
normal and exceptional addresses. Each replay supplies exact values at every nonreserved lookup,
including zero and unique one-bit leaders where a main
correction word is read. An assertion rejects any other
nonreserved access. The exact certificate and fixed-schedule
proof apply to the full tables.
