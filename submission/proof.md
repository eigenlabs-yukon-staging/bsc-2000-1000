# Shared shortened-tail decoder with exact low-bank allocation

The code carries two 256-bit Hamming payload groups, fourteen
32-bit wide groups and a 40-bit tail. Indices start at zero. Channel and machine words are
little-endian. All vector sums below are XOR.

## Code and fixed decoder

For Hamming group g=0 or 1, let its four payload words be
X0, X1, X2, D. Set A=X0 XOR D, B=X1 XOR D, C=X2 XOR D.
Transmit the seven words

    A, B, C, D, A XOR B XOR D, A XOR C XOR D, B XOR C XOR D

at input words 7g through 7g+6. This is an invertible
payload transform followed by 64 binary Hamming(7,4)
lanes. The syndrome columns are (3,5,6,7,1,2,4).

Start with two 16-bit first rows, read bit zero first,
given in construction.py. Pack every cyclic left rotation of
the first row below the corresponding rotation of the second,
for rotations 0 through 15. For columns 16 through 31 pack the
rotations of the transposed second and transposed first rows.
Transposition replaces bit i by bit (-i mod 16). The rows come
from code C64,1 in Anev, Harada and Yankov's Table 1.
[New extremal singly even self-dual codes of lengths 64 and 66](https://arxiv.org/abs/1708.05950)

In each column i flip original parity bit 27 exactly when bit i
of the fixed mask 0x2e27aa03 is one. Let P0(a) be the XOR
of the resulting columns indexed by the set bits of a.
The score calculation uses this exact modified generator.

For the transmitted parity coordinates use P(a)=pi(P0(a)),
where pi puts original bit PARITY_ORDER[j] into physical
position j. The literal order in construction.py is a
permutation of 0 through 31. Thus the wide-component code
is reordered only in its parity coordinates. This choice
changes which eight parity checks are dropped in the
shortened tail, whose protected word uses the top 24
physical bits of P(a).

Put fourteen words (a | (P(a) << 32)) into input words 14
through 27, in payload order, for the fourteen 32-bit
payload groups a. The tail has low 32-bit data word a and
high eight raw data bits b. Put its word

    a | (b << 32) | ((P(a) >> 8) << 40)

into input word 28. All later genuine channel coordinates
are zero at encoding, and undergo the ordinary channel.
This fixes all 2000 transmitted bits. The wide and tail
parts are systematic and the Hamming transform is invertible.
Exactly 2^1000 distinct words are encoded.

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

For the wide code, the syndrome of an error with low and high
halves e and f is P(e) XOR f. For the protected tail, omit the
eight raw bits and use (P(e) >> 8) XOR f on its 56 coordinates.
In each system choose the minimum-weight error with that
syndrome, breaking ties by the lexicographically least
increasing tuple of flipped positions. Parity unit columns
ensure all syndromes have leaders. Denote their low 32 data
bits by Cw and Ct, respectively.

Let rho=0x1739e3b1 and K=2^24. Call r a wide-assigned key
when it is at least K, or when 48<=r<K and its minimum wide
leader weight w for syndrome r XOR rho satisfies

    w <= RESTORE_MAX_WIDE[t],

where t is its minimum tail leader weight for syndrome
r XOR (rho >> 8). The literal thresholds in construction.py
are (-1,-1,-1,-1,4,5,6,10) for t from 0 through 7. The exact
tail histogram proves all its weights fall in that range.
A negative threshold assigns no keys of its tail weight to
wide decoding. Every key from 48 through 2^32-1 stores

| Assignment | High 32 bits | Low 32 bits |
| --- | --- | --- |
| Wide | P(r) XOR rho | Cw(r XOR rho) |
| Tail | (P(r) XOR rho) AND (2^32-256) | Ct(r XOR (rho >> 8)) |

The table generator uses the fixed rank order above to test
the threshold and fill each entry. It uses at most 2^32-48
words, and the underlying input/output region remains disjoint.
The decoder itself takes the same fixed instruction sequence
regardless of the assigned role of any key.

For each wide word load Y, mask it to its low 32-bit half
A, load table/address A, XOR the loaded word into Y,
and shift the result right by 32 to get R. Load at R and
XOR that word into the original Y. The first lookup carries
all 32 required parity bits when A is wide-assigned, or when it
is tail-assigned and the low byte of P(A) XOR rho is zero.
On those first keys, R is the true wide noise syndrome XOR rho.
If R is then a wide-assigned key, the low 32 bits of the result
equal ordinary minimum-distance systematic decoding. For an even-indexed
wide word, mask that value to 32 bits. For an odd-indexed
wide word, shift it left by 32 and OR it into the preceding
even result. This packs the fourteen wide outputs into seven
64-bit output words. Exactly eight instructions are used per
wide word before the seven lane-assembly ORs and output stores.

For the tail, execute the same initial load, low-32 mask,
indirect load, and XOR. Shift right by 40 in this case,
yielding a 24-bit address R. When A>=48, every table entry
gives correct high 24 parity bits, so R is its protected
56-bit error syndrome marked by rho>>8. Load at R and
XOR that result into the original Y. When R is a tail-assigned key, the low
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

    W = 2*24 + 14*8 + 7 + 8 + 16 = 191.

The declared exact comparison checks W^q <= 2^p for
the manifest's exponent p/q. Four illustrative inputs
are replayed on the byte-identical organizer counter, with
actual, exact table entries at every addressed location.
One trace exercises a low wide-assigned key as both a parity
key and a correction key. The fixed schedule proves the
worst-case count directly.

## Error probability

Each Hamming lane succeeds exactly when its noise is the
weight-zero or unique single-bit leader, with probability
h=(9^7+7*9^6)/10^7=531441/625000. There are 128 independent lanes.

In the wide syndrome space let a_w count the wide-assigned
correction keys with minimum leader weight w. There is one
certified correct noise pattern per such key. Write
M=sum_w a_w*9^(64-w). Let D be the number of safe first keys:
all wide-assigned keys, plus the tail-assigned keys whose
unmasked parity byte is already zero. For each exact fixed
noise pattern the received systematic half is uniform on
2^32 values, since the data word is uniform. Thus the lower
bound on wide-component success is

    p_w = M*D/(10^64*2^32).

In the tail let t_w count the tail-assigned correction keys
with minimum tail leader of weight w, and N=sum_w t_w*9^(56-w).
The first lookup only needs the high 24 parity bits, which
are correct at all addresses at least 48, independently of
assignment. The conditional uniform first-key argument again
gives a factor (2^32-48)/2^32. The eight raw channel bits
must remain unflipped independently of the protected noise.
Therefore the tail success lower bound is

    p_t = N*(2^32-48)/(10^56*2^32) * (9/10)^8.

These are sufficient local success events. Their composition
returns the entire payload. Uniform message words and the
disjoint channel noise coordinates make the component events
independent, giving overall success at least h^128*p_w^14*p_t.
Each run checks one minus this exact rational product against
the submitted bound.

## Exact computation of syndrome counts

A mask of weight w signifies w distinct flipped column
positions. The reached set R_w consists of syndromes with
a representation by at most w positions. In the 32-bit
case the declared checker enumerates in the original
parity order, then applies the fixed binary basis map T
in verify.py. The syndrome of a reordered error is pi of
its original syndrome under the corresponding coordinate
permutation. Because this reordering preserves weight,
it maps minimum leader weights exactly. Thus the set of
original syndrome columns suffices for this calculation. Its output bit j is the parity of
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

For the joint allocation counts, first enumerate the shortened
reached sets R_t over its entire 24-bit syndrome space using
the same recurrence. Store the weight t at each newly reached
position. The declared computation reaches all 2^24 positions
and proves this array contains their exact minimum weights.

The wide basis deliberately uses bits 22 through 29 of T(s)
for bits 24 through 31 of the physical syndrome pi(s). This
is checked against the eight explicit basis rows. Consequently
keys r<K, for which the physical syndrome is r XOR rho,
occupy exactly 512 rows of each sector. Their row indices
are ((rho >> 24) << 9) OR j, for 0<=j<512. The checker
constructs T^-1 by eliminating paired transformed and
original unit vectors. For each such row and each of its
2^13 low indices, compute the corresponding original
syndrome by T^-1, apply pi, and XOR rho to get r. The
precomputed linear contributions of the low 13 bits are
reused across the rows. This enumerates all K marked low
keys exactly once across the four sectors.

In each of those rows, among addresses r>=48, partition the
bits into masks for their exact tail weights, obtained from
the shortened array at r XOR (rho >> 8). Store a second
mask within each class for which the low byte of P(r) XOR rho
is zero. The byte contributions to P(r) are computed by
binary linearity on the three bytes of r. At each wide
weight w, the newly reached bits in that row are its weight-w
syndromes. The intersection with each mask gives, by
population count, the joint counts at all (w,t), and their
zero-byte counts. Counting all new bits in the designated
rows also gives the complete marked-low wide histogram,
including the 48 reserved input and output keys.

The computation subtracts that marked-low histogram from
the full wide histogram to get the high-bank correction
mass. For each joint class satisfying w<=RESTORE_MAX_WIDE[t],
add all its keys to the wide correction mass and the
first-key count D, and remove their exact contributions
from the eligible tail leader mass N. For the other classes
add only their zero-byte counts to D. The initial count D
is 2^32-K, and initial N uses all shortened weights minus the
48 direct tail exclusions. All three sums count precisely
the declared sufficient events. The four runs verify the
complete wide counts and joint allocation data separately by
sector. The zero-sector run also checks the fixed schedule
and replays the four illustrative inputs.
