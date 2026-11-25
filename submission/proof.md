# Affine tag decoding with 2 entries per side

## Encoder

Number coordinates starting with the least significant bit. The
first nineteen columns in construction.py define a linear map P
from nineteen bits to thirty-one. The last twelve define a linear map J
from twelve bits to thirty-one, with 4096 distinct images as checked by
the declared experiment. A side with private word u and tag t
is the fifty-bit word E(u,t)=(u,P(u) xor J(t)).

Let G be the binary linear map whose twelve column images
are fixed as GAMMA_COLUMNS in construction.py. The
checker verifies G^2(s)=s xor G(s) and that G(s) differs
from 0 and from s for every nonzero s. For a pair of
twelve-bit parameters (a,s) define the four tags

    (a, a xor s, a xor G(s), a xor s xor G(s)).

For each successive 100-bit message component use the first
four nineteen-bit fields as its private words, the next twelve bits
as a, and the final twelve bits as s. Send E(u,t) for each of its
four sides in tag order. The first nineteen bits of each side
give its private word. Subtracting P of that word from the
parity field and using injectivity of J determines its tag.
The first two tags recover a and s. Thus the encoder is linear
and injective and the ten independent 200-bit components use
all 1000 message bits and exactly 2000 channel bits.

## Decoder and tables

For one received side (U,V), let s=V xor P(U) be its
31-bit private syndrome. The private code of words E(u,0)
is the kernel of this syndrome map. For every syndrome s
choose a minimum-weight binary error l(s), with the
deterministic choice defined by leaders() in construction.py.
Breadth-first search starts from syndrome zero and visits
edges in the fixed coordinate order. Edges XOR a syndrome
column. All syndromes are reachable using parity unit
columns. A shortest path cannot repeat a coordinate,
because two equal flips can be removed. Conversely an
error of weight w gives a path of length w. The breadth-
first records are therefore minimum-weight errors.

Given a hypothesized tag t, the closest side with that
tag is obtained by removing l(s xor J(t)). Its distance
is the weight of that leader. The preloaded leader word
holds this weight in its low eight bits and the private
nineteen correction bits above it.

Extend the twelve independent columns of J to a basis of
the thirty-one-dimensional syndrome space, using parity
units in increasing order and retaining each independent
unit. The echelon routine records the original basis
coordinates of each pivot and split_syndrome returns the
coordinates of a syndrome. Write s=s0(q) xor J(z), with
z in the low twelve coordinate bits and q in the remaining
nineteen bits. Here s0 is the linear map from those
nineteen coordinates to the corresponding basis vectors.
Each y has a unique q and z.

For every q determine two records by sorting
all relative tags d in increasing order of the pair

    (weight(l(s0(q) xor J(d))), d).

One record holds d in its low twelve bits and that weight
above. This finite prescription fixes each value even on
unusual syndromes. The actual side list contains the
absolute tags z xor d with their corresponding preloaded
distances. Indeed J is linear, and

    s xor J(z xor d) = s0(q) xor J(d).

The table generator applies this fixed linear coordinate
split to determine the absolute records for each syndrome.
The online decoder obtains them from the syndrome's
leader word, as described in the counted schedule.

For each pair of distinct sides i,j and their two selected
tags ti,tj, the affine line relation determines all four
tags. To see its explicit form put D=ti xor tj. The
remaining tags are ti xor G(D) and tj xor G(D). The order
of these two tags among the remaining side slots is a
literal function pair_slots(i,j) of the two known slots.
The polynomial identity G^2=G+I gives G^3=I.
The four tag positions are evaluations at 0,I,G,G^2.
Their nonzero differences are I, G and G^2, all invertible.
For the two remaining positions, their normalized
differences from position i are G and G^2. The second
satisfies G^2(D)=G(D) xor D. This proves the formula
also when D=0, when all tags coincide. The checker verifies
each literal pair mapping for all values of the direction.

The decoder tries all six side pairs and all 2
by 2 choices from those two lists. A tag
tuple with its tags present on any two sides is therefore
among the searched candidates. For each candidate it sums
the two distances already stored with the list entries
and the two other leader distances. It keeps the tag
tuple of smallest sum and the first in the fixed loop
order on a tie. There is always at least one hypothesis.
The private output words are the observed private bits
corrected by the appropriate chosen leaders for that
tuple. The tag parameters are a=t0 and s=t0 xor t1, so
the decoder always outputs a valid component message.

## Exact block error certificate

For a received word y and the decoder's output message,
let its distance be the distance from y to that message's
encoded word. Under the uniform message prior, y's
contribution to correct decoding is 9^(200-d)/(2^100*10^200)
if its decoder distance is d. The z-shift calculation
above shows that the shortlist test for a true tag
depends only on the side noise. For a selected quadruple
of private leaders admitted on at least two sides, every
translation by an encoded word therefore gives a
received word with some searched candidate at that
quadruple's weight. The decoder distance is at most
that weight everywhere in the corresponding global coset,
which has 2^100 words.

Let D(t) be a nondecreasing lower bound on the number
of such occupied cosets by radius t. Then there are at
least 2^100*D(t) received words with decoder distance at
most t. By summation by parts with the decreasing
weights 9^(200-t), this gives a lower success numerator

    sum_(t=0..T) (D(t)-D(t-1)) * 9^(200-t),

with D(-1)=0. The omitted contributions beyond T are
nonnegative.

For any finite set of selected errors, let V be its size and
m(s) its multiplicity in a coset. If the number of ordered distinct
pairs within a common coset is at most A0, the number of occupied
cosets is at least the ceiling of V^2/(V+A0) by Cauchy-Schwarz.
For any integer k>=1 a second bound is the ceiling of
(2kV-A0)/(k(k+1)). Indeed for each nonnegative integer m, the
inequality

    k(k+1) * 1[m>0] >= 2k*m-m*(m-1)

is equality for k and k+1 on the positive branch and otherwise
follows from (m-k)(m-k-1)>=0. Summation gives the bound. Both
expressions increase with V for fixed A0 and fixed k. Thus they
also apply with a certified lower value v<=V. The computation
takes k=1+floor(A0/v), the maximum of both ceilings and the
previous lower count. All the volumes used contain zero and v>=1.

For a fixed nonzero word c of weight d in a length n code,
the number of ordered errors (e,e xor c) both at weight at most t
depends only on d. Choose the a positions of e on the support of
c. On that support the two weights are a and d-a. Off of it they
agree, so one can choose at most t-max(a,d-a) of n-d positions.
The exact number is

    H(n,d,t) = sum_a binom(d,a) *
               sum_(0<=b<=min(n-d,t-max(a,d-a))) binom(n-d,b),

with a limited to 0<=a<=d and max(a,d-a)<=t. It vanishes for
d>2t. Summing over all nonzero codewords by their weight counts
gives the pair upper bound, and is exact for the full ball.

For the private code C0, apply the two occupied-volume
inequalities to the local ball and its exact codeword
spectrum. The resulting cumulative D0(t) lower-bounds
the number of private leaders through weight t. The
next section strengthens the profile at low weights and
constructs a lower volume of the quadruples of private
leaders with true tags admitted on at least two sides.
The collision bound below applies to all those quadruples.

The difference of two errors sharing a global
coset is a codeword. If all the difference tags are zero, each
side pair has the same private syndrome and must agree, since
only one leader per syndrome was chosen. Every ordered distinct
pair therefore has a nonzero difference in the two tag parameters.
There are two classes of these differences. For s=0 and a
nonzero the four tags are all a. For s nonzero the four tags
are distinct and either none or one is zero.

For the differences with four nonzero tag slots, drop the
local leader conditions and count pairs using H(200,d,t)
against the full spectrum of codewords with this kind of
tag. For the differences with one zero tag slot, the two
errors in that slot must be the same leader. If that leader
has weight w, there are at most min(binom(50,w),2^31) choices
for it. The remaining three side words of the codeword
difference form a 150-bit word. For weight d of this word
the possible ordered errors in those sides are bounded by
H(150,d,t-w). Multiply by the number of choices for the
common leader and sum over w from 0 to t. Count the three-sided
words separately for each admissible tag tuple. The zero-tag
slot is already accounted for by the identical leader
condition, so this count has no multiplier of 2^19 from
the possible private words in the omitted slot.

For each local tag x write h_x(z) for the local weight
polynomial of E(u,x) over all 2^19 private words. Marking
the total three-side or four-side weights as powers of z,
the tagged spectrum polynomials are products of the
corresponding h_x. For s=0, a nonzero contributes h_a^4.
For s nonzero the sets {0,s,G(s),s xor G(s)} are additive
four-element subspaces. They partition nonzero directions
into triples. Each one is enumerated once by taking the
least nonzero s in it. For that subspace, the four-tag sets
are its additive cosets, with a chosen as the least
member of the coset. Each unordered coset yields twelve
distinct ordered parameter pairs, from the four initial
a choices and the three nonzero s choices, all with the
same polynomial product. The coset containing zero yields
a three-side product of h_s, h_G(s), h_(s xor G(s)),
with the same multiplicity of twelve. All other cosets
yield the four-side products. Together with the constant
nonzero tags these account for exactly 4096^2-1 pairs.

These polynomial sums are computed by packing integer
coefficients into digits of base 2^78. To bound digits let
M=max(1,max_x h_x(1/2)). The evaluation at 1/2 of any product
of up to four side polynomials is at most M^4. Each sum and
each intermediate polynomial in the computation consists
of at most 4096^2 such products, including integer
multiplicities. Since coefficients are nonnegative, its
coefficient at degree d<=50 is at most 2^d*4096^2*M^4.
The checker evaluates the bound at d=50 with exact
denominator 2^50 for each h_x(1/2) and proves it
is below 2^78. Thus ordinary integer products implement
all relevant convolutions without carries into or between
the retained slots.
The needed degree is at most 2t, and t is capped at 25,
so each product and sum can be reduced above degree 50.
The synthetic fourfold convolution of private leader counts uses exact integer shell coefficients.

Before the outer packed products, the checker verifies that
all nonzero tag spectra below weight seven are zero. It
factors out X^7 from each such polynomial. Products of
four factors can then be truncated through degree 22
and products of three through degree 29. Restore the
powers X^28 and X^21, respectively, when extracting
the spectra. The coefficient and digit bounds above
apply to the corresponding original degrees, so the
shifted calculation has the same exact no-carry property.

## Bounding list removals by exact noise weight

Fix a transmitted side tag t and a noise pattern e. The
received private syndrome is J(t) xor syndrome(e).
If the split of syndrome(e) is (q,z_e), the split of the
received syndrome is (q,z_e xor t). Membership of the
true tag t in its actual list is therefore equivalent to
the relative tag z_e being in the preloaded list at q.
This condition depends only on e. If e is a chosen
private leader of its own syndrome, its weight is the
distance used by the table for that relative tag.

Write R for 2. Among rejected private
leaders of weight w in a quotient, each has at least R
distinct preceding relative tags. Their chosen leaders
have weight at most w. Each such ordered pair differs by
a local codeword E(u,t') with t' nonzero, because their
syndromes lie in the same quotient and are distinct. The
weight spectrum of those local words is the sum of the
supplied tag histograms for t' nonzero. The checker bounds
ordered pairs with first weight w and second at most w
even when both errors range over all patterns. For a
nonzero difference of weight d this count is

    sum_a binom(d,a)*binom(50-d,w-a),

where 0<=a<=d, 0<=w-a<=50-d and d-a+(w-a)<=w.
The constraints say that a positions of the first error
lie on the difference support and the second weight is
at most w. Multiply by the nonzero-tag word counts and
sum over d. Call the resulting bound H_w. The number of
rejected leaders at weight w is at most

    min(binom(50,w), floor(H_w/R)).

There is a second cumulative bound. In a quotient with
m private leaders at weight at most t, exactly min(m,R)
of them are admitted by the sorted list. The remaining
(m-R)_+ satisfy, for integers m>=0,

    (4R-2)*(m-R)_+ <= m*(m-1).

For m<R the left side is zero. For m>=R subtract the
left from the right to get (m-2R)*(m-2R+1), nonnegative
at every integer. Sum over all quotients. As before each
ordered distinct pair has a nonzero tag difference. The
full ball pair bounds for the local nonzero-tag spectrum
give an upper count P_t for weight at most t on both
sides. Thus total rejected leaders through t are at most
floor(P_t/(4R-2)).

Through weight eleven the checker counts exact private
leader weights and ranks. It transforms each of the fifty
syndrome columns to coordinates (q,d) using the split
basis, with tag d in the low twelve bits. Form two binary
checks of q using bit masks 8193 and 65553. Move their
outputs to the highest two quotient positions and drop
the former coordinates 0 and 4. The restrictions of the
two masks to these deleted positions have determinant
one over F2, verified by the checker. The remaining
seventeen bits and the two checks thus give a bijective
linear change of quotient coordinates. It preserves each
complete quotient as a row, so it leaves all rank counts
invariant and partitions the rows into four sectors.

To make the row translations shorter, the verifier makes
one further linear change within each quotient row. It
selects the first independent preserving columns under
ordinary high-pivot elimination on the seventeen inner
quotient bits. The pivot stores the reduced inner word
and its tag component. Since these pivots have full rank,
there is a unique linear map F from inner quotient words
to tag words taking all pivot words to their stored tag
components. Replace the tag of each column by its XOR
with F on the inner quotient part. For every fixed
quotient this is an XOR permutation of tag coordinates,
so it keeps all row populations and ranks. In particular
the stored pivot relations cause at least seventeen
preserving columns to have zero tag shift.

A row is an integer bitset, with its bit d marking one
reachable syndrome. Through weight five, enumerate all
increasing coordinate tuples and mark their syndromes in
the designated sector. At weight w from six to eleven,
select min(36,26+3*(w-6)) of the 36 columns that preserve the
sector, ordered by tag bit count with a stable tie order.
Include the preceding rows and their XOR translates by
each selected column, and enumerate all weight-w tuples
from the complementary columns directly. Any weight-w
error using a selected coordinate is covered by the
translation of the error with that coordinate removed,
which is in the same sector. The complement includes
the other weight-w errors. Conversely every translation
can be realized with weight at most w. By induction, the
set in each sector is exact at each of the twelve weights.

For a preserving column XOR its seventeen-bit quotient
part into the row index. Permute tag bit positions by
XOR with its tag part. For tag bit i put step=2^i and
use the low mask (2^4096-1)/(2^step+1). Its left shift
by step selects the complementary lanes. Shift the two
masked pieces in opposite directions by step and OR
them to exchange the index bit. Compose the toggles of
the column's tag part. Each translation reads the
preceding exact rows. Once 7/8 of the sector is reached,
the checker equivalently collects translations into each
destination row, using XOR's own inverse. It may stop
adding terms for a row once all its bits are set.

The population of each row gives its cumulative count
of distinct syndromes. If successive counts are a,b,
the new syndromes have the current weight. Their chosen
minimum leaders occupy ranks a through b-1, independently
of the order among ties. The checker counts rank positions
below 64 with interval endpoint differences and sums the
remaining positions into bucket 64. Each sector checks
its part of ranks.json by this procedure. The four sectors
partition all quotients, so their sum gives all exact
counts by weight and rank through weight eleven. A shortlist
of size R admits precisely ranks below R.

Let T_w lower-bound the cumulative count of private
leaders through weight w and U_w upper-bound the
cumulative count of rejected private leaders. Through the
exact cutoff use the precise rank counts. After it set

    T_w=max(T_(w-1),D0(w)),
    U_w=min(U_(w-1)+min(binom(50,w),floor(H_w/R)),
            floor(P_w/(4R-2))).

The upper recurrence is valid by induction using both
previous inequalities. Its upper bound starts from the
exact earlier rejected count. The actual cumulative
number of admitted leaders is at least T_w-U_w. Let
F_w be the running maximum of this lower bound through w,
also including zero.

Take a synthetic side with T_w-T_(w-1) patterns of each
weight w. Set a_(-1)=0 and put

    a_w=min(a_(w-1)+T_w-T_(w-1), F_w).

Mark a_w-a_(w-1) of these patterns as admitted. The increment
is nonnegative because F_w is nondecreasing. It is at most
the number of new patterns at that weight, so this is an
ordinary finite marked universe. Its cumulative total T_w
and cumulative admitted count a_w are both lower bounds
for the corresponding actual leader counts.

For each fixed choice on the other three sides, consider
replacing the synthetic side by the actual private leaders.
Through a total-weight cutoff, the possible completions on
that side are all leaders within the remaining weight if
the other sides have at least two admitted tags. They are
only admitted leaders if the other sides have exactly one
admitted tag. They are empty if none of the other tags is
admitted. The two cumulative lower bounds imply that each
replacement cannot decrease the number of such quadruples.
Apply this argument to all four sides in succession.

Write A(X) for the synthetic marked shell polynomial and
U(X) for its unmarked shell polynomial. The number of
synthetic quadruples with weight at most t and at least two
admitted tags is the cumulative coefficient through t of

    A(X)^4 + 4*A(X)^3*U(X) + 6*A(X)^2*U(X)^2.

This is a lower bound on admitted volume and contains the
zero error. Use the unchanged upper collision count to
obtain the occupied cosets for each cutoff. Weighted
summation by parts yields a lower component success
numerator N. The ten independently decoded components
have independent uniform inputs and channel flips, so
success is at least N^10/10^2000. The exact comparison
gives block error at most 998973310975542631/1000000000000000000.

## Counted word-RAM schedule

The program is the tuple sequence emitted by schedule() in
construction.py. A string operand denotes a register and an
integer operand denotes a literal. A move is a free copy.
Each other tuple denotes one named bsc-ops-v1 primitive.
All loop bounds, branches and indices of the generator
depend only on the stated code and fixed list size 2.
There is no runtime branch or loop counter.

Addresses 0 through 31 hold the input and 32 through 47
hold the final output. The leader bank begins at 2^31,
contains 2^31 words, and has exclusive end 2^32.
The private parity table begins at 2^30. Its entries
are 2^31 OR P(U), so XOR with the observed parity gives
the marked private syndrome. The following tag table
stores J(t). The pair banks described below are disjoint
from these and the leader table. There are no working
memory arrays.

Every leader entry for syndrome s stores its selected
minimum weight in bits 0..7 and its 19-bit private
correction in bits 45..63. In between it stores the
actual absolute shortlist records for that same
syndrome. Compute its quotient q and tag shift z
offline by the split basis above. For each listed
relative tag d of q, the entry contains d xor z in the
record's low 12 bits and the exact minimum weight in
its top six bits. There are 2 records of 18 bits,
starting at bit 8. Their fields occupy at most 36 bits.
The XOR with z changes only the tag bits. All distances
are at most 31, fitting both listed and leader fields.
This lookup is exactly the same list and correction
mapping as the separate mathematical definitions; the
storage shares one word between the two mappings.
The split basis and quotient records define the table
offline and do not require decoder memory.

Extract each side by loading and shifting its fixed input
word, adding a second load and shifted piece on crossings.
Mask to fifty bits and split into the observed private
and parity portions. Retrieve the marked private syndrome
using the private parity table. Load the leader-bank
entry directly at that address. Shift and mask each
side's selected record to its literal tag and distance
registers. Subsequent leader lookups at adjusted marked
syndrome addresses extract only their distance or
their private correction. The packed records do not
affect those uses.

For each of six side pairs, and every literal pair of
shortlist indices, Read the two known tags and weights as register operands.
There are six read-only banks beginning at 2^27, each
of stride 2^26, with the exclusive end at 2^29. The
banks are below the parity table and above the prefix.
In each bank, use index (t_i<<12) OR t_j in each plane.
Each plane has 2^24 words. Plane 0 gives the pair of tag displacements J(t_k), J(t_l)
packed into two 32-bit lanes. Plane 2 gives all four
completed tags in four 12-bit fields in side order.
Planes 1 and 3 are zero. Table generation uses the
affine identities proved above. The bank base has no bits
below 26 and the two-part key uses only bits 0 through 23.
Form the row prefix once per first tag by shifting and
ORing in the bank base. OR in the second tag to form
the plane-0 key.

For each side pair, pack the two marked syndromes into
the same two lanes once. XOR the pair table word with
that combined value to obtain both marked leader
addresses. Mask and shift to extract the two addresses.
Each leader word stores its distance in the low eight
bits and the private correction above those bits. Add
both whole leader words to the sum of the listed weights.
The four weights are at most 31 each, so their sum is at
most 124 and the final mask to eight bits gives their
exact sum. The higher shortlist and correction fields are multiples
of 256 and do not affect the retained sum. This still
holds if a machine word addition wraps modulo 2^64. The first candidate initializes
the incumbent with its distance and plane-0 address.
Later candidates update these using one unsigned comparison
and two selects. The strict comparison preserves the first
tie. Recover all four tags by reading plane 2 at the winning
address with bit 25 set. The table layout changes none of
the candidates, their scores or their order. After splitting
the four tags, compute t0 xor t1 for the output direction.

Obtain the winning private corrections from the four
leader entries and XOR them into the observed private
words. The last two fields are t0 and t0 xor t1. The
six fields of every component appear in increasing output
bit order. A word beginning at a field boundary starts
with a free move of that field. Each later piece is
shifted and ORed into the partial word. A crossing field
finishes the old word, which is stored, and its shifted
high piece initializes the next word. A field ending on
a word boundary also completes that word. After the final
field the last word is stored as well. Thus each output
word is written exactly once. The computed fields
fit their widths, so unused high bits remain zero.
These output stores cannot affect later input reads or
table access.

The generator first gives every assignment a distinct
value identity. Its operands refer to the most recent
identities for their names at that program point.
This forward pass records the last position using each
identity. The second pass translates the operands using
the current register allocation. It releases a register
only after its final operands for that instruction have
been recorded. An instruction reads its old operands
before assigning a result, so this also permits the
result to reuse a register whose last value it consumes.
Choose the least currently free register for the new
result, and release it after the instruction if it has no
uses. By induction each operand sees the correct earlier
value and different live values have different registers.
This is a generation-time renaming and emits no additional
operations. The verifier checks the finite generated
program uses 33 registers, within the limit
of 64.

Counting the fixed instruction stream gives
760 add, 760 and_, 960 load, 230 ltu, 558 or_, 460 select, 268 shl, 562 shr, 16 store, 370 xor, for exactly W=4944 primitives.
The checked comparison W^100000 <= 2^1227147
establishes L=1227147/100000. The same computation verifies
the affine pair maps, syndrome coordinates and memory
extents. It does not preload the full tables.

## Declared computations

The four supplied local spectrum arrays encode the 4096
tag histograms in groups of 1024. Each uses a Base64-
wrapped zlib stream of unsigned little-endian three-byte
coefficients. The four runs independently enumerate the
2^19 candidate private words for the tags in their sectors
[0,1100), [1100,2200), [2200,3300), and [3300,4096).
They verify all 51 coefficients per tag. In parallel, each
run checks its designated quotient sector of ranks.json
using a separate standard-library process. Both parts
must complete successfully. The fourth run
uses the verified arrays to compute the outer spectra,
pair bounds, shortlist removal volumes and exact weighted
error comparison. It also checks memory boundaries,
paired tag identities and the exact operation comparison.
