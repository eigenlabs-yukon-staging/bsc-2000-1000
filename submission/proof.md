# Eight-sided polynomial tag decoding

## Explicit encoder and side leaders

The encoder splits the message into five 200-bit components.
Each has eight consecutive 19-bit private words followed
by four consecutive 12-bit tag parameters. Throughout,
positions are zero based and the low bit comes first.
All additions of binary vectors and field elements are XOR.

For each side print a private word u in 19 bits and then
P(u) xor J(t) in 31 bits. The 31 literal columns in
construction.py define P by its first 19 columns and
J by its remaining 12 columns. Both maps are binary linear.
The checker verifies that the images J(t) are distinct.
This gives eight sides of 50 bits per component and 2000
bits in the five components.

The tags are evaluations of a polynomial of degree at
most three over the binary field F of order 4096. Use the
polynomial basis with modulus 0x1053, which denotes
X^12+X^6+X^4+X+1. The checker verifies irreducibility by
checking X^(2^12)=X and gcd with X^(2^(12/p))-X for
p in {2,3}, all modulo the stated polynomial as applicable.
The evaluation points are all elements of F satisfying
x^8=x, sorted by their 12-bit binary representations.
The checker verifies there are exactly eight and that
they contain 0,1 and are closed under addition and
multiplication. They form a subfield since every nonzero
element of a finite closed set containing one has its
inverse there. The four tag parameters specify the
values at the first four points. A degree-at-most-three
polynomial with those values is unique by the root bound,
and its evaluations are obtained by Lagrange interpolation.
The eight private words appear explicitly in the output,
and the four parameters are recoverable from the first
four tags given the private words, by injectivity of J.
Thus the encoder is an explicit binary linear injection.

The side syndrome of a 50-bit pattern (u,v) is P(u) xor v.
For each side syndrome s fix the least-weight pattern
l(s), breaking ties by the first discovery in a breadth-first
traversal of the syndrome graph, starting at zero and
trying the 50 coordinate columns in increasing bit order.
An edge corresponds to one XOR flip. The shortest path
length to a syndrome is its minimum pattern weight.
All syndromes occur since the 31 parity unit columns
span them, so the leaders are all determined and have
weight at most 31.

Extend the tag columns to a basis by taking independent
parity unit columns in increasing order. Write each
syndrome as s0(q) xor J(z), where q is a 19-bit quotient
coordinate and z is a 12-bit tag coordinate. For each
q, sort the 4096 relative tags d by

    (wt(l(s0(q) xor J(d))), d).

The first 3 tags form the fixed list for that row.
For a received side with quotient q and shift z, add
z to each listed tag d to get the absolute tag t.
Linear J gives s xor J(t)=s0(q) xor J(d), so the list
records give the exact minimum distances for those tags.
The tie rule is fixed and does not change any equal-weight
shell count.

For each four-element subset of the eight side slots in
lexicographic order, run through the Cartesian product
of the four chosen side lists in rank order. Interpolate
the polynomial from those four selected tags, evaluate
all eight tags and score the sum of the eight local
leader weights. Choose the first candidate of minimum
score. Correct the eight received private words with
the selected side leaders and output the first four
polynomial values as the tag parameters. The lists are
nonempty, so there is always a candidate. For any fixed
tag tuple this rule finds its nearest codeword in the
private coordinates by independent side minimization.
Every tuple with its tag present in at least four
side lists is included in the search.

## Exact local spectra and ranks

For each tag x define h_x[d] as the number of private
words u for which (u,P(u) xor J(x)) has weight d.
These are the full 51-coefficient polynomials for all
4096 tags. Each run computes its designated portion
from the literal columns. It enumerates the 2^19
(u,P(u)) pairs by XOR addition of each private basis
column, starting with zero. For each tag it adds
popcount(u)+popcount(P(u) xor J(x)) once per pair.
The exact coefficient arrays are checked against the
four compressed files. They count all private words
for every tag, including zero.

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

Each quotient sector has exactly 2^29 syndromes. The
checked total of the twelve shell counts in each
sector is 2^29, so all local syndromes and leaders are
accounted for, with no later weights. With R=3
let A_w=sum of the rank counts below R at weight w
and U_w=sum of the rest. These are respectively the
exact admitted and unadmitted chosen leaders on a side.
The weight of each representative is its graph distance.
Its rank among different tags in one quotient determines
whether it is in that quotient's prescribed list. The
order among equal-weight leaders does not affect these
counts. Let A(X), U(X) be their shell polynomials,
and let T(X)=A(X)+U(X). The verifier checks the
admitted total R*2^19.

## Admitted volume and two-error overlap

One component is a binary linear [400,200] code C.
Consider error patterns whose eight side patterns are
chosen private leaders and whose ranks fall in the
shortlist on at least four sides. The set of these
patterns of weight at most t has exact cardinality V_t
given by cumulative coefficients of

    sum_(i=4..8) binom(8,i) A(X)^i U(X)^(8-i).

If such an error occurred from a true codeword, its
rank test on each side is exactly the received shortlist
test for the true tag. Indeed the received syndrome is
the XOR of the true tag displacement J(t) and the error
syndrome. Translation in the tag coordinate preserves
its quotient and takes the relative list to the actual
list. As the errors under consideration are chosen
leaders, the decoder can reproduce the true codeword
among its hypotheses. In particular the decoder finds
distance at most t for all received words in the
corresponding code coset. The assertion is independent
of the true codeword.

An ordered pair of distinct patterns in a given coset
has a nonzero difference in C. On any side with zero
tag difference, the private syndromes of the two local
chosen leaders agree, so those leaders are equal. Thus
that side has zero pattern difference. If the tag
difference is identically zero, all eight patterns are
equal, which is excluded. The tag difference otherwise
comes from a nonzero degree-at-most-three polynomial;
by the root bound it is nonzero on at least five sides.
Let j be the number of nonzero tag positions.

### Cumulative pair counts on one side

Fix a nonzero tag difference x on a changed side.
If the first error is the chosen private leader of
syndrome s, the second leader has syndrome s xor J(x).
This rule gives a bijection, since the same XOR
recovers s. There are T_u leaders of weight u, so
exactly T_u pairs have first weight u before imposing
any condition on the second weight. There are likewise
T_v pairs whose second leader has weight v.

For a fixed difference of weight d, a pair of errors
with first and second weights u,v has
a=(d+u-v)/2 ones on the support in the first error,
and b=(u+v-d)/2 common ones off it. The second
error is complementary on the support. This gives

    binom(d,a) * binom(50-d,b)

ordered pairs, taking zero for nonintegral or invalid
indices. Sum over d weighted by the exact spectrum
h_x[d] of side differences with tag x. Call the result
D_x(u,v). It bounds all candidate leader pairs at
these weights since it also counts other error patterns.

Let C(w)=sum_(u<=w) T_u, with C(w)=0 for w<0.
For each first weight u and combined-weight cutoff w,
the second weight is at most w-u. The cap for each v
is min(T_v,D_x(u,v)). Their sum over eligible v can
also be capped by T_u, by uniqueness of the second
leader from the first. Therefore a cumulative upper
bound on ordered pairs of total weight at most w is

    F_x(w)=sum_u min(T_u, sum_(v<=w-u) min(T_v,D_x(u,v))).

There is another upper bound for every integer k.
If the total is at most w, the first weight is
at most k or the second is at most w-k-1. Otherwise
their sum would be at least w+1. The union gives
C(k)+C(w-k-1). Take the minimum across k from -1
through the largest local leader weight, and with F_x(w).
Each bound is nondecreasing in w. Define K_x[w]
as consecutive increments of this minimum, starting from
zero at cutoff -1. K_x gives a nonnegative synthetic
shell sequence with cumulative counts at least the
actual local pair counts.

At a fixed tuple of tag differences, replace each
changed side's pair measure by its K sequence. This
increases or preserves the combined pair count through
any cutoff. To check one replacement, fix the weights
on all other sides and apply the cumulative dominance
to the remaining weight. Sum over the other weights
with their nonnegative multiplicities, and repeat for
each side. We may therefore use K as shell coefficients
for the rest of the support calculation.

### Five changed sides

For a fixed support of five, a nonzero degree-three
polynomial that vanishes on the other three points is
a nonzero scalar times the monic product of those three
linear factors. The checker verifies that the affine
maps x -> a*x+b of the eight-point subfield, with a
nonzero, take the first three sorted field points onto
every triple of points. Such maps permute the eight
positions. Substituting an affine map in a polynomial
takes evaluation vectors to their permuted vectors.
It also permutes the nonzero multiples of the monic
vanishing polynomial. All sides use the same local
pair sequences, so every triple of zeros gives the same
combined-weight polynomial after summing over scalars.

For that representative triple, let r_i be the five
values of its monic vanishing polynomial on the changed
positions. Let K_x(X) have coefficients K_x[w]. Sum
the exact synthetic products

    sum_(a != 0) product_(i=1..5) K_(a*r_i)(X).

Each unchanged side gives the same error twice, so its
combined-weight polynomial is T(X^2). Multiply the
sum by T(X^2)^3 and compute its cumulative coefficients
through 2t. They upper-bound the pairs with both weights
at most t for this fixed support. Since the local caps
are at most 2^31 in total (take threshold k=-1), the
total sum of product coefficients is at most
4095*(2^31)^8. Base 2^288 therefore gives exact
carry-free packed arithmetic through degree 108.

### Six, seven or eight changed sides

Fix a changed support S of size j. The tag polynomial
vanishes at its 8-j fixed complementary points. It
factors as their monic vanishing polynomial times a
polynomial of degree less than m=j-4. Each m-position
projection on S is injective on the possible tag tuples:
the vanishing factor is nonzero at those positions, so
polynomial interpolation determines the quotient.
Each changed position has a nonzero tag difference.
For fixed totals w_i of the two local error weights, the
synthetic product for a fixed tuple of tags uses the
factor K_(t_i)[w_i] for each changed side. The cumulative
comparison just proved justifies using these coefficients
in an upper count after the combined cutoff. We omit the
shortlist admission conditions when counting these pairs.

For a fixed ordered tuple of totals w_i define
g_i(x)=K_x[w_i] on nonzero x. For each m-element subset
I of S set H_I=product_(i in I) g_i(t_i)^(j/m).
Each index occurs in m/j of these subsets, so the
geometric mean of the H_I is product_i g_i(t_i).
Finite Hölder's inequality and injectivity give

    sum_(tag tuples) product_i g_i(t_i)
      <= product_I (sum_(tag tuples) H_I)^(1/binom(j,m))
      <= product_i (sum_(x != 0) g_i(x)^(j/m))^(m/j).

For the first inequality, when all sums are nonzero,
normalize each H_I by its sum, use the arithmetic-
geometric mean inequality at each tuple, and sum over
tuples. If a sum is zero the left side is zero. The
second inequality follows by extending each subset
sum to the entire Cartesian product of nonzero tags,
since the projection is injective.

Use fixed scale S0=2^24. For each w compute

    I_w=sum_(x != 0) ceil_root_m(S0^m * K_x[w]^j),
    Q_w=ceil_root_j(S0^(j-m) * I_w^m).

Here ceil_root_k gives the least nonnegative integer
whose kth power reaches the argument, by exact binary
search. The first expression bounds the inner norm sum
after multiplying it by S0, and Q_w/S0 is therefore
an upper bound on the complete norm. Define
Q(X)=sum_w Q_w*X^w. The scaled polynomial Q(X)^j
bounds the ordered changed-side products by combined
weight, after dividing by S0^j.

On each unchanged side the two chosen leaders are
identical, so their total weight has polynomial T(X^2).
Use the full product Q(X)^j*T(X^2)^(8-j), divide by
S0^j and round upward after summing its coefficients
through degree 2t. For each pair in which both
component error weights are at most t, their combined
weight is at most 2t. Thus the result is an upper count
for that support, even though it allows a larger
set at the combined cutoff. Multiplying by binom(8,j)
and summing j=5 through 8 gives an upper count P_t
for all distinct admitted pairs in the same coset.

Each product is computed with base 2^512, with the
total at X=1 checked to be smaller than that base.
Thus the products and the cumulative sums use exact
coefficients and can discard terms above degree 108.
The checked file profiles.json contains the cumulative
integer counts for a fixed support of each size.
Each support profile is verified in its designated
run before the final score uses the complete array.
The admitted-volume polynomial also fits base 2^512
because the total number of local leader products is
at most (2^31)^8.

### From pairs to component success

Let P_t be the resulting upper total count of ordered
distinct retained pairs in the same component coset.
If a coset contains n retained patterns, its occupancy
indicator satisfies, for integer k>=1,

    k*(k+1)*1[n>0] >= 2*k*n-n*(n-1).

At n=0 both sides vanish; at positive n subtract
to get (n-k)*(n-k-1) >= 0. Sum over cosets and
bound using V_t and P_t. The checker uses
k=1+floor(P_t/V_t) and rounds the quotient up.
It also takes the lower bound
ceil(V_t^2/(V_t+P_t)) by Cauchy-Schwarz, as verified
by writing V_t as the sum of positive occupancies.
The running maximum of both gives a valid integer
lower bound D_t on the number of successful
distance-at-most-t cosets, since the set of received words with
decoder distance at most t grows with t.
Here V_t never vanishes since the zero error is
admitted on every side.

All received words in an occupied coset have some
admitted candidate at distance at most t. Over uniform
200-bit messages and for each received word, correct
decoding has probability contribution
9^(400-d)/(2^200*10^400), where d is the distance
to the codeword output by the decoder. Each coset has
2^200 words. Using D_(-1)=0 and summation by parts
against decreasing shell probabilities, the component
success is at least N/10^400, where

    N=sum_(t=0..54) (D_t-D_(t-1))*9^(400-t).

The five disjoint message components are independent
under the uniform message and their channel flips are
independent. The program below implements this component
rule separately on each of them. Hence the block error is
at most 1-N^5/10^2000. The declared computation verifies
its exact upper rational score 99169299528205687/100000000000000000.

## Counted 64-bit RAM schedule

The literal tuples emitted by schedule() specify the
decoder. A string is a named register and an integer
is a literal. A move is a free copy. Each other entry
is exactly the named primitive. All decisions in code
generation depend only on fixed parameters or literal
indices. There is no runtime branching or implicit
address, loop, call or return work.

The preloaded coefficient rows start at 2^25. There
are 70 groups and four rows per group, one for each
chosen side. A row has four planes of 4096 words.
The offset within a plane is one 12-bit tag. From
the Lagrange basis polynomial for this selected side,
compute its contribution at every evaluation point.
Planes 0 and 1 store J of the four complementary
contributions in two pairs of 32-bit lanes, one pair
per word. Both halves have their top bit clear.
Planes 2 and 3 hold the tag contributions themselves,
four 12-bit fields in each word. XOR across the four
selected sides therefore gives exactly the required
polynomial evaluation, or its displacements by linear J.
With the row base aligned to four planes, OR with the
tag and the plane bits forms the indexed address.

At 2^26 the descriptor table has one four-word record
for each group and each literal rank tuple. The four
words contain, in respective low and high halves,
the working list address for that side/rank and the
coefficient row address for its plane 2. Every such
word is determined by the literal group and rank indices.
The finite ranges checked by the verifier put both banks
above word 47 and below the parity table, disjoint.

The parity table begins at 2^30 and stores P(u) with
the bit for 2^31 set. The tag table stores J on every
tag value. Lower and upper coordinate tables then
store the linear coordinate transformation of all
16-bit low and 15-bit high segments of side syndromes.
The quotient bank contains 3 entries per quotient,
each a relative tag and its exact weight above 12.
A small writable array has one absolute list of that
size per side. The read-only leader bank begins at 2^31
and extends through word 2^32-1. Its word for each
syndrome contains the leader weight in the low byte
and its 19-bit private correction above. This defines
all preloaded values in the finite RAM. All ranges
are disjoint and checked.

Decode a component by loading and stitching words as
necessary to isolate each successive 50-bit side. Mask
out its private bits and shift to obtain parity bits.
XOR the marked parity table word into the observed
parity to get the marked private syndrome. Use its
two parts for the split tables. The high part already
contains the marker bit at position 15, accounted for
by subtracting 2^15 in its literal table base.
XOR their outputs for the quotient and tag shift.
Load each prescribed relative record, XOR the shift
into its tag field, and store it in the fixed working
list. No candidate depends on a preceding component's
working lists, since all its own entries are written
before they are read.

For each chosen group pack the four complementary
marked syndromes into two 64-bit words, in the same
lane order as planes 0 and 1. The generator expands
every rank tuple in the stated order. The fixed tuple (0, 3, 3, 3) gives the number of list
entries cached for each of the four chosen slots in order.
The cached entries are those with the lowest literal ranks.
Each cached entry supplies its weight and two packed
contributions. Other entries are looked up as needed.

At the first node XOR the two contribution words into
the packed marked syndromes. The next three nodes
XOR the two contribution words for the next selected
tag into their running packed words, sharing prefix
work. The known-side weights are added in parallel to
this accumulation, with the first such weight copied
to initialize its sum. Each leaf then already holds
the marked leader addresses for the four complementary
sides. No runtime tree traversal or counters enter the
emitted schedule.

Mask and shift each pair to extract the two addresses.
Load both leader words and add them in full to the four
selected weights.
The high correction fields are multiples of 256.
All eight local weights are at most 31, so their
sum is at most 248. Mask to one byte to get exactly
that sum, irrespective of the correction fields.
The first candidate initializes the best distance
and descriptor address. Later candidates use an
unsigned strict comparison and two selects to update
them. Ties therefore keep the first candidate.

Recover the winning eight tags by the four descriptor
words at the winning literal address. For each selected
side use its stored tag to read planes 2 and 3 of the
specified row. XOR both words into two 48-bit accumulators.
The first is the output tag parameter block and both
supply the eight winning absolute tags. For each side
look up J of its tag and XOR with the marked syndrome
for its selected leader address. Shift the leader word
right by eight and XOR the private correction into
the received private data. This implements the same
nearest-private reconstruction used in the error proof.

Pack the eight private corrections followed by the
four tag parameters into successive 19- and 12-bit
fields in the output words. Each shift and OR is an
emitted primitive. A word beginning at a field boundary
starts with a free copy, and a crossing field writes
its low part into the completed word and shifts its
remaining high bits into the next word. For a component
boundary inside an output word, store the partial word
and reload it at the start of the next component.
Those stores and loads are included in the count.
The final partial word has only 40 significant bits
since the fields have total length 1000. All upper
24 bits are zero.

Each expanded component instruction stream is assigned
physical registers at generation time. Give each virtual
destination a new value number and track the last
instruction reading each numbered value. A result
may reuse an operand's register after the operand's
last read, since a primitive reads all its operands
before writing its output. Take the least currently
free register for each live result. Induction over
the stream gives the same operands and results with
distinct simultaneous values in distinct registers.
Allocate each component separately since its only carried
value is an explicitly stored and loaded output word.
The verifier checks the emitted stream uses 61
registers, within 64, and counts exactly

    154610 add, 89510 and_, 126474 load, 28345 ltu, 40 mul, 9268 or_, 56690 select, 788 shl, 61177 shr, 140 store, 84320 xor.

Thus W=611362. The exact comparison
W^100000 <= 2^1922167 proves L=1922167/100000.

## Declared computations

The four runs independently check all 4096 tag spectra
in contiguous groups [0,1100), [1100,2200), [2200,3400),
[3400,4096), and verify their own four disjoint quotient
sectors of the rank data. Each run executes the two
independent tasks in parallel using standard-library
processes and asserts they both completed. Each spectrum
file is a Base64-wrapped zlib stream of unsigned
little-endian three-byte counts. The last run also
checks the exact outer occupancy, the error ceiling,
the field modulus and interpolation maps, the
memory extents, and the fixed instruction count and
exponent. The check does not materialize the complete
leader or coefficient banks.
