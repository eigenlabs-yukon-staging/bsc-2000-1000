# Affine tag decoding with 4 entries per side

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

For every q preload 4 records by sorting
all relative tags d in increasing order of the pair

    (weight(l(s0(q) xor J(d))), d).

One record holds d in its low twelve bits and that weight
above. This finite prescription fixes each value even on
unusual syndromes. The actual side list contains the
absolute tags z xor d with their corresponding preloaded
distances. Indeed J is linear, and

    s xor J(z xor d) = s0(q) xor J(d).

Splitting a syndrome is a linear operation. The low-split
table contains the coordinates of each 16-bit syndrome
suffix. The high-split table contains the coordinates of
each 15-bit word at positions 16 through 30. XOR the
two outputs to obtain the coordinates of the full syndrome.

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

The decoder tries all six side pairs and all 4
by 4 choices from those two lists. A tag
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

## Exact shortlist admission volumes

Fix a transmitted side tag t and a noise pattern e. The
received private syndrome is J(t) xor syndrome(e).
If the split of syndrome(e) is (q,z_e), the split of the
received syndrome is (q,z_e xor t). Membership of the
true tag t in its actual list is therefore equivalent to
the relative tag z_e being in the preloaded list at q.
This condition depends only on e. If e is a chosen
private leader of its own syndrome, its weight is the
distance used by the table for that relative tag.

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

In each weight shell w the checked table counts the
private leaders by rank. All four sectors reach all
2^31 syndromes in total by weight eleven. Thus there
are no further weights to bound. Let T_w be the exact
number and A_(i,w) the exact count below the fixed list
size of side i. Ranks tied at a weight still have exact
prefix counts because the list order resolves every tie.
Set U_(i,w)=T_w-A_(i,w).

Sum the eleven products of four shell polynomials whose
admitted side sets have size at least two. For a side in
the set use A_i and for any other use U_i. These products
are disjoint by their exact admitted side set. Cumulative
coefficients give the retained volume V_t, without any
unknown multiplicities. The count includes the zero error.

## Pair bounds on occupied cosets

Two distinct retained leader quadruples in one coset
have a codeword difference with a nonzero tag tuple.
If all tags agreed then each local leader would be
the same because a syndrome determines one chosen
leader. Any two affine tag coordinates determine both
parameters, so at most one difference tag is zero.
The possible changed support sizes are three and four.

For a fixed tag difference x, the exact local spectrum
h_x[d] counts its codeword differences by Hamming weight.
For a word difference of weight d and ordered error
weights u,v, the first error has a=(d+u-v)/2 ones on the
difference support and b=(u+v-d)/2 elsewhere. The ordered
ambient count in that cell is

    D_x(u,v)=sum_d h_x[d] binom(d,a) binom(50-d,b),

where invalid integer choices contribute zero. This also
bounds the count with both errors restricted to be chosen
leaders. Adding J(x) bijects their syndromes. Consequently
those leader pairs use any leader on either side at most
once and the respective shell marginals are at most T_u
and T_v. On an unchanged side the pair must be the same
leader, whose shell is exactly T_u.

For total weight cutoff w, one upper bound counts all
pairs with first weight at most k or second weight at
most w-k-1. Its value is C(k)+C(w-k-1), where C is the
cumulative T count. Another upper bound charges all second
weights at most k once as C(k), then uses for each first
weight u the smaller of T_u and the sum of min(T_v,D_x(u,v))
over k<v<=w-u. Minimize both bounds over fixed k values
-1 through eleven. All terms before minimization are
nondecreasing in w. The cumulative upper bound therefore
has nonnegative increments K_x[w]. The verification loop
can stop increasing k when C(k) reaches the current cap,
since all residual terms are nonnegative. For fixed tag
differences the convolution of these synthetic local
sequences cumulatively dominates the true total weight
sum. This follows by successively replacing each local
sequence and summing against nonnegative coefficients.
The unchanged side uses T_u at total pair weight 2u.
A pair with each error weight at most t has sum at most 2t.

A separate bound keeps both weight coordinates. For each
fixed first weight u, the cumulative count through second
weight v is bounded by

    H_(x,u)(v)=min(T_u, sum_(z<=v) min(T_z,D_x(u,z))).

Take its consecutive nonnegative increments as a row in
a polynomial in X,Y, with X recording u and Y recording
the second weight. At any fixed vector of first weights,
convolving the rows preserves cumulative dominance for
the second weight. Summing over all first-weight vectors
of total at most t thus gives an upper count for the two
cutoffs. For an unchanged position use T_u at X^u Y^u.

### Summing over the affine tags

Fix a support of size j. For m=j-2, projection of its
possible tag tuples onto any m chosen changed positions
is injective. The zero positions together with those m
positions give two known coordinates, which determine the
tag tuple. For fixed local exponents let f_i(x) be the
nonnegative synthetic coefficient on changed side i.
For each m-subset E of changed side indices, take the
product of f_i(x_i)^(j/m) over E. The geometric mean over
these subsets equals the full product over changed sides.
Hölder bounds its sum over the tags by the geometric mean
of the sums for each E. Injectivity bounds each such sum
by the product of the individual sums of f_i(x)^(j/m).
Counting how many subsets include each position gives
the per-coordinate norm (sum_(x!=0) f_i(x)^(j/m))^(m/j).

For the sum-cutoff calculation use this norm on each
total-weight coefficient. At scale S=2^24 use the integer

    ceil_root_j(S^(j-m) * (sum_x ceil_root_m(S^m f_i(x)^j))^m).

The nested upward roots make this at least S times the
norm. The j changed polynomials use these coefficients
and the unchanged polynomial uses exact common weights.
Sum coefficients through degree 2t and divide by S^j
rounding up. The integer product uses digits of base
2^256. Its whole coefficient sum, the product value at
1, is bounded by the digit mask in the verifier, so no
coefficient carries can occur.

For support four the same norm and upward roots apply
to each pair (u,v) instead of a sum of weights, using
scale S=256. Sum coefficients having both exponents at
most t and divide by S^j rounding up. The product uses
carry-free integer packing with exponent position
u*45+v. There is no exponent aliasing because the
largest possible summed second weight is 4*11=44.
The coefficient base uses enough bytes to hold the
product evaluation at (1,1). Truncating after first
weight 25 preserves all coefficients needed.

For support three there is a sharper exact sum of
synthetic products. Fix the single zero side. For each
nonzero direction s, the zero coordinate determines a.
The remaining tags are exactly the three members of
{s,G(s),s xor G(s)} in some order. The checked identity
G(G(s))=s xor G(s), with neither zero nor fixed nonzero
images, partitions the 4095 directions into 1365
triples. Each triple contributes the same product for
its three nonzero directions. Sum those products of
the three two-weight polynomials and multiply by three
and the common leader polynomial. Both exponent cutoffs
then give the support-three bound for any fixed zero side.
The same exponent stride 45 applies. The verifier bounds
the base digits by three times the number of triples times
the largest changed coefficient sum cubed times the common
coefficient sum. This bounds the total of all coefficients,
so all additions and products have exact digits.

The witness stores the minimum of the summed-weight and
the separate-cutoff counts for each t and each support size.
The overall pair cap is four times the support-three
profile plus the full support-four profile.

Use the resulting pair count P_t and the exact V_t in the
occupancy inequalities above. The weighted summation by
parts at T=25 gives lower success numerator N for one
component. The ten independently decoded components have
independent uniform data and channel flips. Success is at
least N^10/10^2000. The exact comparison gives block error
at most 199340877440341041/200000000000000000.

## Counted word-RAM schedule

The program is the tuple sequence emitted by schedule() in
construction.py. A string operand denotes a register and an
integer operand denotes a literal. A move is a free copy.
Each other tuple denotes one named bsc-ops-v1 primitive.
All loop bounds, branches and indices of the generator
depend only on the stated code and fixed list size 4.
There is no runtime branch or loop counter.

Addresses 0 through 31 hold the input and 32 through 47
hold the final output. The leader bank begins at 2^31,
contains 2^31 words, and has exclusive end 2^32.
The private parity table begins at 2^30. Its entries
are 2^31 OR P(U), so XOR with the observed parity gives
the marked private syndrome. The tag table stores J(t),
the gamma table stores G(t), and the split tables and
quotient list bank follow. The list bank has 2^19 rows
of 4 words. Every shortlist tag and distance remains in its own named register.
All used table and working ranges are disjoint and lie
above the input/output prefix. Table values are fully
determined by table_word.

Extract each side by loading and shifting its fixed input
word, adding a second load and shifted piece on crossings.
Mask to fifty bits and split into the observed private and
parity portions. Retrieve the marked private syndrome
using the private parity table. The upper portion includes the marker as bit 15. Subtract 2^15 from the high table base in its address literal to account for that bit. XOR the two
split table entries to obtain the twelve-bit tag shift
and nineteen-bit quotient. At the literal list indices,
load the selected relative tags and weights and convert
the tags with the tag shift. The fixed constructor emits
only the stores required by the tag placement above.

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
exact sum. The high correction fields are multiples of
256 and do not affect the retained sum. The first candidate initializes
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
program uses 49 registers, within the limit
of 64.

Counting the fixed instruction stream gives
3160 add, 2280 and_, 3320 load, 950 ltu, 40 mul, 1398 or_, 1900 select, 388 shl, 1362 shr, 16 store, 1290 xor, for exactly W=16104 primitives.
The checked comparison W^50000 <= 2^698757
establishes L=698757/50000. The same computation verifies
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
must complete successfully. Sectors 1 and 2 also independently verify the support-three
and support-four cumulative profiles from profiles.json.
Sector 3 combines these bounds with the exact admission
volumes and weighted occupancy, and verifies the error
and work inequalities, memory and pair identities.
