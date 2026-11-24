# Word-aligned 64-bit components in 454 operations

## Code and local syndrome leaders

Use zero-based binary indices, low bit first. Divide
the first 992 message bits into 31 consecutive blocks
of 32 bits each. For each block a, send a followed
by its 32-bit parity P(a). To define P, take the
binary polynomial G specified by the integer
141459731, which has degree 27. For each data
bit i=0..31 form the polynomial word with integer
representation

    u_i = (1 << (27+i)) | (EXTRA[i] << 59).

The four-bit EXTRA values are explicit in
construction.py. Its parity column consists of the
27 low remainder bits of u_i modulo G, followed
by the four bits EXTRA[i], then the parity of all
the bits of u_i plus its remainder. The terms of
u_i and its remainder have disjoint positions.
XOR the columns for the nonzero bits of a to
obtain P(a). This fixes a linear systematic
64-bit component of dimension 32.

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

## Word RAM and reserved addresses

Let rho=OFFSET and H=(2^32)-3*2^16. The preloaded
correction word at address r in 48..H-1 is
C(r XOR rho), with high bits zero. After this range
are three consecutive 65536-word tables. At H+j
store the parity of the low sixteen data bits j,
with the fixed rho XORed in. At H+65536+j store
the parity of the high sixteen data bits j. At
H+2*65536+j store the chosen eight-bit tail decision
for received block j, shifted left by 32. This is
the table_word rule. These four ranges are disjoint
and use only the available addresses below 2^32.
Addresses 0..31 contain the input, and addresses
32..47 remain zero until the final stores.

The machine program is the literal expansion of
schedule(). Per input word j=0..30 load the word,
mask its data part and shift to obtain V. Isolate
the low and high sixteen data bits using a mask
and a shift. OR the low value with H and the high
value with H+65536 to obtain the two parity
lookup addresses. Their low sixteen base bits are
zero. Load the two words and XOR both into V.
This computes the correction address r=s XOR rho.
Load at r and XOR into A. Thus all thirteen
listed scalar primitives have fixed execution and
compute the ordinary result when r is safe.
The two parity addresses and r are always within
the memory bound.

For an even j move the result directly to
output accumulator U[j//2]. For an odd j shift
it left by 32 and OR into that accumulator.
This uses two additional primitives on each of
the fifteen odd indices. On safe keys each
result occupies exactly its intended 32-bit lane.
On an exceptional key the lookup may read input
or zero prefix words or a word from one of the
three top tables. The resulting bits can affect
other data positions, but the actual program is
defined for every such input.

The main program finishes with U15's low lane
assigned by j=30. Load word 31, which contains
only the sixteen received tail bits, OR in
TAIL_BASE, load its table word and OR it into U15.
Mask U15 to its low forty bits. This guarantees
all 24 padding bits are zero even on exceptional
inputs. Store all sixteen words to words 32..47.
The count is 31*13+15*2+4+1+16=454, with the
comparison 454^1000<=2^8827 checked exactly.
There are sixteen persistent accumulators and
six reusable registers. All branches and choices
in schedule() concern literal indices during code
generation. The emitted program contains no runtime
loop, branch, call or implicit memory operation.

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
sectors. Fifty-one transformed columns have top
two bits zero and thus preserve all four sectors
under XOR. The other thirteen are available
for direct enumeration but are never used in
a sector-local bitset translation.

For a sector let R_t be its transformed syndromes
attained by errors of weight at most t. Enumerate
all masks through weight five, keeping only
syndromes in that sector. For each weight t>=6
choose a subset E of the preserving columns
using the smallest low-bit counts in the exact
sort order of the checker. The counts selected
are 33,36,40,44,51 at weights 6,7,8,9 and all
later weights, respectively. Begin with R_(t-1),
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
entire syndrome space and all 48+3*2^16
excluded keys, without assuming an
unverified histogram. Each run checks
saturation at weight twelve.

Let A be the sum of 9^(64-w) over the
safe syndrome leaders counted in these
four histograms. Let T be the sum of
9^(16-w) over the tail leaders, which the
checker enumerates completely by BFS.
For every message, safe main success has
probability A/10^64 and tail success has
probability T/10^16. All these local events
depend on disjoint independent channel
flips. On their intersection the counted
program packs the right data in every lane.
Thus block error is at most
1-A^31*T/10^2000. The exact rational
comparison gives the claimed
999960630397147047/1000000000000000000.
Sector zero's experiment also checks the
instruction count and uses the byte-identical
reference counter on sample received words
including all categories of exceptional key.
The all-input cost and output format
follow from the fixed expansion above.
