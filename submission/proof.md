# Ten-sided polynomial tag decoding

## Explicit encoder and side leaders

The encoder splits the message into four 250-bit components.
Each has ten consecutive 19-bit private words followed
by five consecutive 12-bit tag parameters. Throughout,
positions are zero based and the low bit comes first.
All additions of binary vectors and field elements are XOR.

For each side print a private word u in 19 bits and then
P(u) xor J(t) in 31 bits. The 31 literal columns in
construction.py define P by its first 19 columns and
J by its remaining 12 columns. Both maps are binary linear.
The checker verifies that the images J(t) are distinct.
This gives ten sides of 50 bits per component and 2000
bits in the four components.

The tags are evaluations of a polynomial of degree at
most four over the binary field F of order 4096. Use the
polynomial basis with modulus 0x1053, which denotes
X^12+X^6+X^4+X+1. The checker verifies irreducibility by
checking X^(2^12)=X and gcd with X^(2^(12/p))-X for
p in {2,3}, all modulo the stated polynomial as applicable.
The evaluation points are the distinct field elements
represented by the integers 0 through 9. The five tag
parameters specify the values at the first five points.
A degree-at-most-four polynomial with those values is
unique by the root bound and its evaluations are obtained
by Lagrange interpolation. The ten private words
appear explicitly in the output and the five parameters
are recoverable from the first five tags given the private
words, by injectivity of J.
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

For each five-element subset of the ten side slots in
lexicographic order, run through the Cartesian product
of the five chosen side lists in rank order. Interpolate
the polynomial from those five selected tags, evaluate
all ten tags and score the sum of the ten local
leader weights. Choose the first candidate of minimum
score. Correct the ten received private words with
the selected side leaders and output the first five
polynomial values as the tag parameters. The lists are
nonempty, so there is always a candidate. For any fixed
tag tuple this rule finds its nearest codeword in the
private coordinates by independent side minimization.
Every tuple with its tag present in at least five
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

One component is a binary linear [500,250] code C.
Consider error patterns whose ten side patterns are
chosen private leaders and whose ranks fall in the
shortlist on at least five sides. The set of these
patterns of weight at most t has exact cardinality V_t
given by cumulative coefficients of

    sum_(i=5..10) binom(10,i) A(X)^i U(X)^(10-i).

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
difference is identically zero, all ten patterns are
equal, which is excluded. The tag difference otherwise
comes from a nonzero degree-at-most-four polynomial;
by the root bound it is nonzero on at least six sides.
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

Further cap this cumulative count by separating a prefix
of second-error weights. Fix k between -1 and the largest
local leader weight. All pairs with second weight at most
k use distinct second leaders and contribute at most
C(k). For the remaining pairs with first weight u and
total weight at most w, sum the caps min(T_v,D_x(u,v))
over k<v<=w-u. Clip each first-weight sum at T_u using
uniqueness of its second leader. This gives the upper
bound

    C(k) + sum_u min(T_u,
                     sum_(k<v<=w-u) min(T_v,D_x(u,v))).

The k=-1 case is F_x(w). Another bound for each k is
C(k)+C(w-k-1), since otherwise the sum of the two
weights exceeds w. Take the minimum of both bounds
over this one fixed set of k values. Each bound is
nondecreasing in w, so the minimum is nondecreasing.
The implementation starts from the best two-threshold
bound and stops scanning k once C(k) reaches the
current cap, since the remaining terms are nonnegative
and C(k) is nondecreasing in k. The result is unchanged.
Define K_x[w] as consecutive increments of that final
minimum, starting at zero before weight zero. This is a
nonnegative synthetic shell sequence whose cumulative
counts dominate the corresponding actual local pairs.

At a fixed tuple of tag differences, replace each
changed side's pair measure by its K sequence. This
increases or preserves the combined pair count through
any cutoff. To check one replacement, fix the weights
on all other sides and apply the cumulative dominance
to the remaining weight. Sum over the other weights
with their nonnegative multiplicities, and repeat for
each side. We may therefore use K as shell coefficients
for the rest of the support calculation.

### Combining the changed support

Fix a changed support S of size j. The tag polynomial
vanishes at its 10-j fixed complementary points. It
factors as their monic vanishing polynomial times a
polynomial of degree less than m=j-5. Each m-position
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
Use the full product Q(X)^j*T(X^2)^(10-j), divide by
S0^j and round upward after summing its coefficients
through degree 2t. For each pair in which both
component error weights are at most t, their combined
weight is at most 2t. Thus the result is an upper count
for that support, even though it allows a larger
set at the combined cutoff. Multiplying by binom(10,j)
and summing j=6 through 10 gives an upper count P_t
for all distinct admitted pairs in the same coset.

Each product is computed with base 2^640, with the
total at X=1 checked to be smaller than that base.
Thus the products and the cumulative sums use exact
coefficients and can discard terms above degree 134.
The checked file profiles.json contains the cumulative
integer counts for a fixed support of each size.
Each support profile is verified in its designated
run before the final score uses the complete array.
The admitted-volume polynomial also fits base 2^640
because the total number of local leader products is
at most (2^31)^10.

### Bounding each of the two error weights

A pair in the cutoff set has each component error weight
at most t. The preceding bounds instead use their sum
at most 2t. We get another upper bound by tracking
both weights. For each nonzero tag difference x,
fixed first leader weight u and second cutoff v, put

    H_(x,u)(v) = min(T_u, sum_(z<=v) min(T_z,D_x(u,z))).

This upper bounds cumulative pairs in that row of the
weight array. Each second leader is matched to a unique
first leader, so each individual shell pair is at most
T_z as well as D_x(u,z). Each first leader is used once,
so the sum is capped at T_u. The nonnegative increments
from H_(x,u)(-1)=0 give a synthetic row whose cumulative
counts dominate the actual row for that fixed first
weight. The program obtains the increments by visiting
z in increasing order and filling the remaining T_u
capacity by at most the shell bound at each step.

For a fixed tag difference tuple and fixed first-weight
tuple, convolve the synthetic rows in the second weight.
Cumulative dominance through any second cutoff holds by
successively replacing each actual row, since convolution
with nonnegative counts preserves cumulative dominance.
Then sum over first-weight tuples of total at most t.
On each unchanged side the corresponding polynomial
is exact, with coefficient T_u at X^u Y^u. This yields
a valid upper bound for the two separate cutoffs.

For j changed sides and interpolation freedom m=j-5,
apply the same projection and Hölder argument above
to each separately specified pair of weights (u,z).
On a changed side use scale 256 instead of 2^24 for
the nested upward roots. The resulting bivariate norm
polynomial has denominator 256 per factor. Raise it
to j and multiply by the unchanged polynomial to
power 10-j. Sum all coefficients with both exponents
at most t, then divide by 256^j rounding upward.
This bound is used for supports of size six through nine.
The sum-cutoff bound alone is used for support size ten.

For the integer product encode X^u Y^z at position
u*(10*11+1)+z. The maximum possible second exponent
is 10*11, so there is no exponent aliasing between
rows. Use enough whole bytes per digit to hold the
total of all coefficients of the product, computed
from the polynomial values at (1,1). This prevents
digit carry. Reduction at the end of the row for
first exponent 67 preserves all needed coefficients.
The computation takes the lesser of the two integer
cumulative bounds at every cutoff for each of these support
sizes, before combining the supports and applying
the occupied-syndrome argument.

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
250-bit messages and for each received word, correct
decoding has probability contribution
9^(500-d)/(2^250*10^500), where d is the distance
to the codeword output by the decoder. Each coset has
2^250 words. Using D_(-1)=0 and summation by parts
against decreasing shell probabilities, the component
success is at least N/10^500, where

    N=sum_(t=0..67) (D_t-D_(t-1))*9^(500-t).

The four disjoint message components are independent
under the uniform message and their channel flips are
independent. The program below implements this component
rule separately on each of them. Hence the block error is
at most 1-N^4/10^2000. The declared computation verifies
its exact upper rational score 49064807638756941/50000000000000000.

## Counted scalar RAM schedule

The literal tuples emitted by schedule() have a named
register for each string and a literal for each integer.
A move is a free copy. Every other tuple gives exactly
one named primitive, with divu and modu denoting
unsigned division and remainder. All control flow in the generator
concerns literal groups, ranks, offsets, and word operands.
The generated stream has no implicit branches or calls.

Two read-only field banks start at 2^24 and 2^25, each
of 2^24 words. A subrow for field coefficient c contains
J(c*t) in the first bank and c*t in the second at the
index t. The lower twelve bits of the subrow start are
zero, so ORing in t implements the address exactly.

The group bank at 2^26 gives a finite mapping for each
of the 252 five-element subsets in lexicographic order.
A group has five rows, each of length eleven. The first
word in a row gives the start of its side's working list.
The other ten give the start of the corresponding
coefficient subrows of the second field bank. They are
determined by the explicit Lagrange formula, one address
for each of the ten evaluation points. All bank ranges
are disjoint and checked.

The marked private parity bank begins at 2^30 and stores
P(u) OR 2^31. The tag table stores J(t). The next two
tables give contributions to the coordinate split from
the low sixteen and high fifteen bits of the syndrome.
The quotient bank holds 3 relative list records per
quotient, as twelve-bit tags followed by exact weights.
The working arrays store the received private words,
marked syndromes, absolute lists, outgoing contribution
words and weights, and two words for the incumbent.
The read-only leader bank begins at 2^31. Each word
stores its private leader weight in the low byte and the
chosen nineteen-bit private correction above. Its
exclusive end is 2^32. The layout is explicit in the
generator and all used ranges are checked.

Setup stitches words to isolate each of the ten 50-bit
sides. Store its private word and its marked syndrome,
using the marked parity table. Use the split tables
to obtain q,z. The upper portion includes the marker
bit, whose value 2^15 is subtracted from that bank's
literal address base. Read the prescribed quotient list
row, XOR z into the tag field of each entry, and store
the absolute records in the side's list.

For each selected group, pack the first four missing
syndromes in two 64-bit pairs. Put the fifth syndrome
in the low half of a third word with a zero high half.
This costs 4+4+1 operations. For each selected side/rank
record, load its tag and weight and store the weight
in a four-operation sequence. Use its Lagrange
coefficients to read J of its field contributions.
Each complete pair needs two address ORs and loads, a
shift by 32, an OR and a store, totaling seven operations.
The single contribution only needs an address OR, load
and store. Thus group setup costs
9+5*R*(4+7*2+3) operations if all entries use
working storage. The last 3 literal selected positions
keep each rank's weight and three contributions in named
registers. Those positions skip four stores per rank.
This saves 4*3*R operations in group setup.

The fixed rank tree expands all R^5 candidate paths.
At depth zero, load its next local weight and three
contribution words, copy the weight and XOR the
contributions into the three marked words. At later
depths add the new weight to its parent prefix and
XOR each contribution into its respective parent word.
Thus at depth zero each node costs seven operations,
and at each further depth it costs eight when fetched
from working storage. At a cached depth the four
operands are already in registers, saving four read
operations per node. The number
of nodes at depth d is R^(d+1). The identical rank
prefix is used for every descendant in that subtree.

Each leaf has the five missing side leader addresses.
For the two complete pairs the mask, shift, two loads and
two additions cost six operations each. The fifth address
uses one load and one addition.
A leader word has the exact local weight in its low byte
and the private correction in the higher bits. Add each
unmasked word to the known-side weight sum. The checked
quotient rank spectra cover all local syndromes by weight 11,
so the total weight is at most 10*11=110<256.
All additions preserve the weight sum modulo 256. One
final byte mask therefore yields the exact distance.
The first candidate initializes the
incumbent by free copies. Later candidates cost a
strict comparison and two selects. A tie leaves the
earlier candidate. Thus ordinary leaves cost 18
operations and the first is three cheaper.

The stored incumbent key is g*R^5+r, where g is its
group and r is the base-R code of its rank tuple.
The group ends with two stores of the key and smallest
distance. Each group after the first has two loads of
these values before its candidates. Consequently the
first group costs 6158, and the remaining groups
cost 6163 each, using the explicit formula above.

For the winner read the group by integer division and
read its literal base-R digits by division and remainder.
Use the group bank row pointers to read the tags for
its chosen side/rank records. Read every tag contribution
from the second field bank, then XOR the five
contributions at each evaluation point. This recovers
all ten tags and outputs the first five as information.
For each side compute the leader address from its tag
and marked syndrome, read the selected private correction
and XOR it into the received private word.

Pack ten nineteen-bit private words and five twelve-bit
tags in order for each component. The free copy rule for
aligned starts and explicit shifts and ORs for crossings
put each field in its designated output bits. Store a
partial output word when a component ends inside it,
and load that word when the next component begins.
These actions are counted. The final partial word
has forty data bits and twenty-four upper zeros.

The register allocator numbers each virtual result and
finds its last use within one stage. It reads operands
before assigning the result and releases registers at
last use, preserving all operands by induction. Each
setup, group, or reconstruction stage receives its own
allocation, because values carried across stages are
explicit memory data. The group shape and the named
register lifetimes depend on R and the fixed packing,
with a separate first group shape. Changing the group
or its side indices changes only address literals and
other constants. The verifier checks representative
schedules for both group shapes and generates all four
setup and reconstruction stages. At most 62 registers
are used in any stage.

The setup/reconstruction operation pairs are
(340, 394), (344, 399), (344, 396), (340, 399). Accounting for the first group and
251 ordinary groups per component gives exactly
W=6215240. The checker verifies
W^50000 <= 2^1128369, establishing L=1128369/50000.

## Declared computations

The four runs verify all tag spectra in contiguous
intervals [0,960), [960,1960), [1960,3040), [3040,4096).
The stored coefficient arrays use zlib and Base64 with
unsigned three-byte coefficients. Every run also verifies
its disjoint quotient-rank sector, using a separate
standard-library process so the two calculations can
proceed in parallel. A successful rank process exit
is required. The runs recompute support profiles for
sizes 6, 7, 8, and 9 and 10, respectively.
The last run verifies the precise occupancy and error
bounds, field and table formulas, fixed stage shapes
and work exponent. It does not materialize the entire
preload banks.
