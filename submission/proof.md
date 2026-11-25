# Shared shortened-tail decoding in 147 operations

The code has three 256-bit payload groups, six 32-bit payload
groups, and a 40-bit tail. Indices start at zero. Channel and
machine words are little-endian. All vector sums below are XOR.

## Code and fixed decoder

For one of the three 256-bit groups, call its payload words
X0, X1, X2, D, in that order. Set A=X0 XOR D, B=X1 XOR D,
C=X2 XOR D. Transmit the seven words

    A, B, C, D, A XOR B XOR D, A XOR C XOR D, B XOR C XOR D.

Use input words 0 through 20 for the three groups in order.
The four payload words can be recovered from the first four
transmitted words. Independently in each of the 64 bit
positions, this is the Hamming code with parity-check
columns (3,5,6,7,1,2,4), represented as three-bit integers.

Let P be the 32-to-32-bit linear parity map whose ith
column, for i from 0 to 31, is generated as follows.
The degree-27 binary polynomial G has integer representation
141459731, with coefficient of X^j at bit j. The 27 low bits
of the column are the remainder modulo G of the binary
polynomial represented by (1 << (27+i)) OR (EXTRA[i] << 59).
Its high five bits are EXTRA[i]. The literal tuple EXTRA
in construction.py supplies all 32 values. For a word a,
P(a) is the XOR of the columns indexed by its one bits.

Put six words (a | (P(a) << 32)) into input words 21
through 26, in payload order, for the next six 32-bit
payload groups a. The tail has low 32-bit data word a and
high eight raw data bits b. Put its word

    a | (b << 32) | ((P(a) >> 8) << 40)

into input word 27. All later genuine channel coordinates
are zero at encoding, and undergo the ordinary channel.
This fixes all 2000 transmitted bits. The wide and tail
parts are systematic, and the Hamming payload transform is
invertible. Thus exactly 2^1000 distinct words are encoded.

### Hamming groups

Write S0,S1,S2 for the syndrome bit planes. Their expressions
are

    S0 = (A XOR D) XOR B XOR P0
    S1 = (A XOR D) XOR C XOR P1
    S2 = (B XOR D) XOR C XOR P2

on the received bit planes, where P0,P1,P2 are the last
three received planes of the group. For each lane, the
syndrome equals the syndrome of the noise. Its Hamming
leader is the single-bit error in the coordinate labelled
by the syndrome if nonzero, or zero if the syndrome is zero.
Correcting that leader yields a valid Hamming codeword.
The first payload plane A XOR D changes on single
errors at labels 3 and 7 only, so its correction mask is
S0 AND S1. The analogous masks for B XOR D and C XOR D
are S0 AND S2 and S1 AND S2. For D the mask is
S0 AND S1 AND S2. Word logic evaluates these functions
independently over all lanes. The fixed template in
construction.py takes seven loads, nine XORs to form the
three raw payloads and the syndrome, four ANDs including
the shared pair in the triple mask, and four XOR corrections.
These total 24 ops per group before the output stores.

### Shared table

For the wide code, define the syndrome of a 64-bit error
with low and high halves e and f to be P(e) XOR f. For the
protected part of the tail, omit the raw byte and take the
syndrome of its 56-bit error with data and parity
subwords e,f to be (P(e) >> 8) XOR f. In each syndrome
system choose the least-weight error, ties broken by
lexicographically least increasing tuple of flipped positions.
The parity unit columns span the respective syndrome
spaces, so a leader always exists. Write Cw and Ct for
the low 32 data bits of the corresponding leaders.

Let rho=0xf7c90230 and K=2^24. At each address r from
48 through 2^32-1 preload a word whose high and low
32-bit halves are as follows.

| Address | High half | Low half |
| --- | --- | --- |
| 48 <= r < K | (P(r) XOR rho) AND (2^32-256) | Ct(r XOR (rho >> 8)) |
| K <= r < 2^32 | P(r) XOR rho | Cw(r XOR rho) |

The two banks share a parity field except for its bottom
eight bits in the lower bank, which are zero so that a tail
correction leaves its raw byte unchanged. The table is wholly
specified by construction.py and the fixed leader order.
It uses 2^32-48 RAM words, and the disjoint 0..47 region
contains the input and output. No scratch memory is needed.

For each wide word load Y, mask it to its low 32-bit half
A, load table/address A, XOR the loaded word into Y,
and shift the result right by 32 to get R. Load at R and
XOR that word into the original Y. If A and R are both
at least K, R is the wide noise syndrome XOR rho and the
low 32 bits of this corrected value are ordinary
minimum-distance systematic decoding. For an even-indexed
wide word, mask that value to 32 bits. For an odd-indexed
wide word, shift it left by 32 and OR it into the preceding
even result. This packs the six wide outputs into three
64-bit output words. Exactly eight instructions are used per
wide word before the three lane-assembly ORs and output stores.

For the tail, execute the same initial load, low-32 mask,
indirect load, and XOR. Shift right by 40 in this case,
yielding a 24-bit address R. When A>=48, both table banks
give correct high 24 parity bits, so R is its protected
56-bit error syndrome marked by rho>>8. Load at R and
XOR that result into the original Y. When R>=48 the low
32 bits are corrected using Ct and the raw received byte
is preserved. Mask the result to its low 40 bits. This
takes eight operations before its final store.

All indirect addresses are 32-bit or 24-bit integers on
every input. If an address is below 48, its read is an
ordinary memory load. Raw inputs have been preloaded, and
no output stores have occurred at the time of any read,
so the output-region values are still their initial zeros.
Even on these paths the last mask clears all output
padding. construction.py fixes the complete literal
instruction schedule, with all generating loops and
indices resolved before the scored program executes.
There is no runtime control flow or implicit address
arithmetic. The 16 output registers and at most 23
additional named registers fit in 64 registers. The
constant operation count is

    W = 3*24 + 6*8 + 3 + 8 + 16 = 147.

The declared exact comparison checks W^q <= 2^p for
the manifest's exponent p/q. Three illustrative inputs
are replayed on the byte-identical organizer counter, with
actual, exact table entries at every addressed location.
Their sampled counts do not supply the worst-case argument.

## Error probability

A Hamming lane succeeds when its actual noise is its
selected leader, hence with probability
h=(9^7+7*9^6)/10^7=531441/625000. These 192 lanes
use disjoint independent noise.

In the wide syndrome space let a_w be the number of
syndromes with minimum leader weight w whose correction
keys are at least K. Exactly one noise word per syndrome
is the selected leader. Let M=sum_w a_w*9^(64-w).
This counts certified noise mass with denominator 10^64.
For fixed noise on one wide component, the received
systematic half is the uniform payload XOR a fixed noise
half. Therefore its parity lookup key is uniform on
2^32 addresses conditional on each noise word. The additional
safe key constraint A>=K multiplies its certified noise
probability by (2^32-K)/2^32. The resulting wide component
success lower bound is

    p_w = M*(2^32-K)/(10^64*2^32).

Similarly let t_w be the number of 24-bit shortened
syndromes with minimum leader weight w whose marked
correction addresses are at least 48, and set
N=sum_w t_w*9^(56-w). The tail first lookup condition
A>=48 occurs with probability (2^32-48)/2^32 conditional
on each fixed noise word. The raw eight bits need zero
flips, independently of the protected component. Thus
the tail success lower bound is

    p_t = N*(2^32-48)/(10^56*2^32) * (9/10)^8.

These are sufficient local success events regardless of
the decoder's behavior on discarded addresses. The message
fields are independent under the uniform message, and
their channel coordinates are disjoint. Since the assembly
on these events returns exactly the decoded payloads,
the whole-block success probability is at least h^192*p_w^6*p_t.
The declared experiment uses the exact integers and tests
that one minus this product is at most the manifest's B.

## Exact computation of syndrome counts

A mask of weight w signifies w distinct flipped column
positions. The reached set R_w consists of syndromes with
a representation by at most w positions. In the 32-bit
case the declared checker applies the fixed binary basis
map T in verify.py. Its output bit j is the parity of
the input syndrome AND row mask BASIS[j]. The elimination
test verifies that T has rank 32, and therefore preserves
reachability and minimum weight when applied to every
column. The top two T bits partition the space into
four disjoint sectors. In the shortened case no basis
change or sector split is used.

In a given wide sector, explicitly enumerate through
weight five. For each following weight w, choose the
listed number of the columns that preserve the sector.
These columns have both top bits zero after transformation.
Start with R_(w-1), OR in its translation by each selected
column, and enumerate all weight-w tuples over the
complementary column positions, keeping the ones within
the sector. If an error of exactly weight w uses a
selected position, deleting it leaves a weight-(w-1)
error in the same sector. This shows that the updates
include all reachable syndromes. Conversely a translation
adds at most one flip, so each newly marked syndrome is
reachable. The complement enumeration is exact.
By induction the resulting set is precisely R_w.
The shortened histogram uses the same algorithm on its
whole space, with direct enumeration through weight four
and the remaining position counts given in verify.py.

The bitset representation divides each sector into rows
indexed by high local bits, with packed bits for the
low local index. The wide row width is 2^13 bits and
the shortened row width is 2^10 bits. Translation by
a column XOR permutes rows by its high part. For each
set low bit j of the column, the corresponding low-index
permutation interchanges lanes of length 2^j using
the alternating mask ((2^(2^low))-1)/(2^(2^j)+1)
and its shift by 2^j. The union uses only the preceding
reached set as input, with full rows allowed to skip
further updates. At saturation every syndrome of
the sector, or of the shortened space, is reached.
Differences of consecutive population counts are
minimum-weight counts.

For the wide exclusions, the key R<2^24 means
the top byte of its original syndrome equals rho>>24.
The inverse binary map T^-1 is constructed by eliminating
paired transformed and original unit vectors. The map
q(z)=T^-1(z)>>24 is a linear projection to a byte.
For each low-index bit position the checker constructs
its q value and groups the positions by that value into
256 masks. For a row, whose sector and high local bits
are fixed, the excluded positions are those whose
low-index q value XOR the q value of the fixed bits
is rho>>24. Thus the row's exclusion mask is exactly
one of the stored masks. For the shortened exclusions,
the forty-eight original syndromes r XOR (rho>>8) are
marked directly. Intersect each reached set with its
exclusions before taking consecutive differences.
This computes the excluded histogram alongside the
full histogram, so subtracting gives a_w and t_w.
Each of the four wide histograms is recomputed by
its corresponding admitted run. Sector zero also
recomputes the entire shortened histogram, checks the
Hamming calculation and counts and replays the schedule.
All four runs check the exact product bound using the
histogram literals whose entries they jointly certify.
