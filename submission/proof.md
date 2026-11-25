# Twenty-sided polynomial tag decoding

## Explicit encoder and side leaders

The encoder splits the message into two 500-bit components.
Each has twenty consecutive 19-bit private words followed
by ten consecutive 12-bit tag parameters. Throughout,
positions are zero based and the low bit comes first.
All additions of binary vectors and field elements are XOR.

For each side print a private word u in 19 bits and then
P(u) xor J(t) in 31 bits. The 31 literal columns in
construction.py define P by its first 19 columns and
J by its remaining 12 columns. Both maps are binary linear.
The checker verifies that the images J(t) are distinct.
This gives twenty sides of 50 bits per component and 2000
bits in the two components.

The tags are evaluations of a polynomial of degree at
most nine over the binary field F of order 4096. Use the
polynomial basis with modulus 0x1053, which denotes
X^12+X^6+X^4+X+1. The checker verifies irreducibility by
checking X^(2^12)=X and gcd with X^(2^(12/p))-X for
p in {2,3}, all modulo the stated polynomial as applicable.
The evaluation points are the distinct field elements
represented by the integers 0 through 19. The ten tag
parameters specify the values at the first ten points.
A degree-at-most-nine polynomial with those values is
unique by the root bound and its evaluations are obtained
by Lagrange interpolation. The twenty private words
appear explicitly in the output and the ten parameters
are recoverable from the first ten tags given the private
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

The first 7 tags form the fixed list for that row.
For a received side with quotient q and shift z, add
z to each listed tag d to get the absolute tag t.
Linear J gives s xor J(t)=s0(q) xor J(d), so the list
records give the exact minimum distances for those tags.
The tie rule is fixed and does not change any equal-weight
shell count.

For each ten-element subset of the twenty slots, in
lexicographic group order, use only the rank tuples whose
sum of literal rank costs is at most 6. The costs
for ranks 0 through 6, in order, are (0, 1, 2, 3, 4, 5, 6).
Each such tuple specifies ten tags, whose unique
polynomial of degree at most nine is interpolated.
Evaluate all twenty tags and score the sum of the
twenty chosen private leader distances. Keep the first
minimum-score tuple. Its corrected private fields and
first ten tag values are the decoded message. The
all-rank-zero tuple is searched in every group. For
fixed tags the independent side corrections give a
closest word consistent with those tags.

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
accounted for, with no later weights. These exact shell
counts classify every relative rank used by the code.
Ranks beyond the stored list are never examined.

## Admitted volume and two-error overlap

One component is a binary linear [1000,500] code C.
Consider noise patterns whose twenty local errors are
the chosen leaders for their own syndromes. Each such
leader has an exact quotient-relative rank. Translating
by the transmitted tag sends this rank to the matching
absolute entry without changing its cost. Rank beyond
the list is treated as infinite cost. The true tuple
is searched exactly when the sum of the ten least finite
costs, counting multiplicities, is at most 6.
Indeed the ten smallest costs minimize over every ten-side
interpolation group. This invariant is independent of the
true codeword.

For p a finite cost, let E_p(X) be the shell polynomial
of leaders at that cost and H_p(X) the shell polynomial
at larger costs, including infinite cost. Let L_p(X,Y)
sum the shells at lower costs, with Y-degree the cost.
Distinguish the actual tenth-smallest cost p, which is
unique even with ties. There are exactly j<=9 sides of
strictly lower cost and at least 10-j sides of cost p.
For each such j the lower-cost j sides contribute
L_p(X,Y)^j, restricted to Y-degree at most
6-(10-j)*p, then evaluated at Y=1. If the
cutoff is negative the contribution is zero. Choose
their positions in binom(20,j) ways. In the other
20-j positions choose k>=10-j equal-cost sides in
binom(20-j,k) ways and multiply their E_p(X)^k by
H_p(X)^(20-j-k). This counts every admissible pattern
exactly once by its tenth cost and its actual sets of
lower and equal positions. The total's cumulative weight
coefficients are the admitted volume V_t. They include
the zero error.

The verifier packs X-degree into consecutive 624-bit
digits. To also track Y-degree, place X^w Y^j at index
w+221*j, since 20*11=220 is the largest total weight.
Each coefficient and its total count are at most
(2^31)^20, checked against the digit mask. This
prevents digit carries in any of these products and
exponent aliasing. Truncating the Y degree at the
remaining cap discards only larger costs, and selecting
the resulting digit blocks and adding them gives
the required sum. The final sum uses the original
exact rank witnesses and nonnegative integer counts.

If two distinct retained patterns share a coset their
difference is a codeword. Its tag difference cannot be
all zero since one local syndrome selects one local
leader. A nonzero degree-at-most-nine polynomial has
at most nine roots, so the changed side support has
size at least eleven. For a fixed support the omitted
zero tags and any j-10 of its j changed tags determine
the whole difference polynomial. The pair bound below
counts all retained pairs, dropping the rank-budget
condition to obtain valid upper counts.

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
vanishes at its 20-j fixed complementary points. It
factors as their monic vanishing polynomial times a
polynomial of degree less than m=j-10. Each m-position
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
Use the full product Q(X)^j*T(X^2)^(20-j), divide by
S0^j and round upward after summing its coefficients
through degree 2t. For each pair in which both
component error weights are at most t, their combined
weight is at most 2t. Thus the result is an upper count
for that support, even though it allows a larger
set at the combined cutoff. Multiplying by binom(20,j)
and summing j=11 through 20 gives an upper count P_t
for all distinct admitted pairs in the same coset.

Each product is computed with base 2^1280, with the
total at X=1 checked to be smaller than that base.
Thus the products and the cumulative sums use exact
coefficients and can discard terms above degree 270.
The checked file profiles.json contains the cumulative
integer counts for a fixed support of each size.
Each support profile is verified in its designated
run before the final score uses the complete array.
The admitted-volume polynomial also fits base 2^1280
because the total number of local leader products is
at most (2^31)^20.

### Combining separate error cutoffs using residues

An admitted pair at distance cutoff t has both individual
component error weights at most t. We bound that event
using a second convolution. For each nonzero tag shift
x and exact first leader weight u, the number of pairs
through second leader cutoff v is at most

    H_(x,u)(v) = min(T_u, sum_(z<=v) min(T_z,D_x(u,z))).

This follows from the two shell marginals and the ambient
shell pair count derived above. The program visits z in
increasing order and fills the row's remaining T_u up to
the cap on each shell. The resulting increments of H
are nonnegative. For each fixed first-weight tuple,
successively replacing actual second-weight rows by these
synthetic rows preserves cumulative dominance, since
convolution with nonnegative counts does so. Sum over
all first-weight tuples with sum at most t. The unchanged
side's exact polynomial has T_u at X^u Y^u. Thus the
synthetic products dominate the true square cutoff.
After this comparison we bound those synthetic products
by coefficient envelopes and then relax their cutoff.

For j changed sides, use the projection and Hölder
argument above coefficientwise at every separate weight
pair (u,z), with m=j-10 and scale S=4096. For each pair let
Q be the resulting integer upper bound on S times its
norm. Before multiplying, replace Q by ceil(Q/2^(u+z)).
An unchanged side's T_u is replaced by
ceil(S*T_u/2^(2*u)). If the total of both weights is d,
multiplying the corresponding product coefficient by
2^d/S^20 restores an upper bound. This is because each
local rounded entry times its own weight power divided
by S bounds the original entry, coefficient by coefficient.

To shorten the polynomial product, encode each monomial
with first index the sum of its two weights and second
index the second weight modulo 16. Sum coefficients
that share these indices. The product therefore retains
the exact total weight d and the residue of the second
weight. A genuine pair with both weights at most t must
have d<=2t and its second weight in the integer interval

    max(0,d-t) <= second weight <= min(t,d).

Include exactly the residues represented by that interval,
each once. Any term whose two weights obey the genuine
cutoff is counted, and additional included terms only
increase the upper bound. The verifier multiplies the
summed allowed coefficients by 2^d for each d and divides
their total by S^20, rounding up at the end.

For exact integer multiplication, put the digit for
coefficient (d,r) at position d+441*r. The largest
total-weight degree of any complete product is
20*(2*11)=440. Thus the degree never wraps into the
residue coordinate. All multiplication carries are
excluded by using enough bytes per digit to hold the
sum of all coefficients of the full product, determined
at (1,1). This also bounds every intermediate product
used in exponentiation because the factors have positive
integer coefficient sums. Folding at the digit boundary
441*16 implements the residue relation modulo 16.
Indeed, the high digits represent second-weight residues
at least 16; adding them into the low positions replaces
the index by its residue. Nonnegative digits and the
same coefficient-sum bound prevent carries after folding.
This yields the exact integer product in the stated
degree and residue representation.

The declared computations take the lesser of this bound
and the original summed-weight upper bound separately at
each cutoff for each support size before applying the
occupancy calculation.

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
500-bit messages and for each received word, correct
decoding has probability contribution
9^(1000-d)/(2^500*10^1000), where d is the distance
to the codeword output by the decoder. Each coset has
2^500 words. Using D_(-1)=0 and summation by parts
against decreasing shell probabilities, the component
success is at least N/10^1000, where

    N=sum_(t=0..135) (D_t-D_(t-1))*9^(1000-t).

The two disjoint message components are independent
under the uniform message and their channel flips are
independent. The program below implements this component
rule separately on each of them. Hence the block error is
at most 1-N^2/10^2000. The declared computation verifies
its exact upper rational score 190084835817278329/200000000000000000.

## Counted scalar RAM schedule

The literal instruction stream is emitted by schedule()
in construction.py. Strings name registers and integers
are literals. A move is a free copy. Every other tuple
is one primitive, with divu and modu denoting unsigned
division and remainder. The generator expands only
fixed loops over literal groups, prefixes, ranks and
bit offsets. There is no runtime loop or branch work
outside the emitted primitives.

The two field lookup banks begin at 2^24 and 2^25
and have length 2^24 words each. At offset c*4096+t
they store J(c*t) and c*t, respectively. The group
bank at 2^26 gives ten rows of twenty-one words for each
ten-element side subset in lexicographic order. In a row
the first word is its selected side's working-list start,
and the remaining twenty words give the aligned subrow
addresses of its Lagrange coefficients in the second
field bank. The explicit field formula fixes all values.

The private parity bank at 2^30 contains 2^31 OR P(u).
It is followed by the tag table of J(t), the two linear
split tables for low sixteen and top fifteen syndrome bits,
and 7 list records per quotient row. A list record
contains a relative tag in the low twelve bits and its
exact leader weight above. The writable arrays store
private observations, marked syndromes, absolute lists,
selected weights and packed group contributions, the
five marked group words, and the incumbent distance
and key. The read-only leader bank starts at 2^31.
A leader word contains its minimum weight in the low byte
and the chosen private correction in the next 19 bits.
The exclusive end is 2^32. Every fixed bank and working
range is disjoint with checked endpoints.

Initialize each component by extracting each 50-bit side
from the input words. Store its received private bits
and form the marked syndrome by the private parity
lookup XOR into its received parity field. Use both
split tables to obtain its quotient and tag shift.
The address for the upper segment subtracts 2^15
because that segment includes the marker. Read the
relative records, XOR the tag shift into each record's
low tag field and store its absolute record.

For a selected group, pack the complementary ten marked
syndromes into five words, each with two 32-bit lanes.
The two loads, shift, OR and store use five operations
per word. For each eligible selected-side/rank entry,
load the entry, separate its tag and weight and store
the weight. For each of five pairs, read J of the two
needed field contributions from the first bank, using
each fixed coefficient and the tag as the key. Two OR
addresses, two loads, a shift, an OR and a store use
seven operations per pair. Group preparation therefore uses
25+10*R*(4+7*5)=2755 operations for R=7.

For each group, enumerate all prefixes of 4 ranks that
meet the total cost cap, in their literal order. The remaining
rank choices of each such fixed prefix are expanded into an
independently allocated scalar stage. At a node, the generator
only emits the branches whose next cost fits the residual cap.
Neither the cap nor this traversal decision is runtime data.

A stage loads the five saved marked pair words, and in an
ordinary stage also the earlier incumbent distance and key. Its
fixed prefix has one node per depth, with eleven primitives
at depth zero and twelve at every later depth. Each free node
loads its own six contribution words, XORs five displacement
words into its prefix and adds the weight. For each given
residual budget, convolve the finite rank-cost list once per
remaining depth. Its cumulative counts through that budget
give the number of emitted nodes at that depth, including
leaf nodes at the final depth. Every descendant gets the same
explicit prefix values. They remain in registers with the
fixed lifetime allocation.

Each fixed-prefix stage preloads the weight and five
displacement words for each of the first 3 ranks at the last
position whose cost fits its remaining budget. Only these literal
ranks can be reused and the stage pays six loads for each of them.
Every leaf whose last rank is cached omits those six loads.
To compute the savings, use the count by summed cost at the
penultimate depth. For a prefix with consumed cost u, count the
cached ranks whose cost is at most the stage budget minus u and
multiply by the number of those prefixes. Subtract six times
this count, add six times the preload count, and leave all
other node and leaf charges unchanged.

At a leaf the ten complementary marked leader addresses
are in five packed words. Split each pair, load the two
leader words and add both, six operations per pair.
A leader word has the exact local weight in its low byte
and the private correction in the higher bits. Add each
unmasked word to the known-side weight sum. The checked
quotient rank spectra cover all local syndromes by weight 11,
so the total weight is at most 20*11=220<256.
All additions preserve the weight sum modulo 256. One
final byte mask therefore yields the exact distance.
The first leaf initializes from free copies.
Later leaves use a strict unsigned comparison and two
selects, preserving the first minimum on a tie. The
stored key is g*R^10+r, with g the group and r its
literal rank tuple in base R. The stage stores its incumbent words before the next stage.
The first stage saves two loads and three first-leaf update
operations. All other stages have counts determined by the
remaining budget and the literal number of prefix nodes.
For each remaining budget, the prover generates a representative
ordinary stage and checks both its primitive count and
liveness. The first stage has its own checked shape. Changing
a prefix that has the same cost only changes address and key
literals, not the virtual register schedule. Multiplying each
stage count by the number of such prefixes in one group gives
472307 operations including group preparation
and before the one five-operation saving per component.

At the end read the winning group and its rank digits
using counted unsigned division and remainder. Use the
row pointers in its group bank to read the ten winning
absolute tags. The second field bank provides each
of their Lagrange contributions. XOR the contributions
to obtain every winning tag. Read each corresponding
private leader using its marked syndrome and J of its
tag, take the correction bits and XOR into the private
observations. Pack the twenty corrected nineteen-bit
words and the first ten twelve-bit tags in order. The
generator emits the shifts, masks and ORs along every
word crossing. The first component stores its final
partial word and the second loads it. The last word
has only 40 data bits and the upper 24 bits are zero.

The generator allocates registers independently for each
preparation, fixed-prefix stage, and reconstruction.
All values carried between them are in the explicit
working arrays. Assign each virtual result its own value
number and record its last use. Read operands first
and return their registers after the last use. The
result is assigned a free register. Induction on
the generated instructions preserves every virtual value
used by an operand. The preparation shape is the same in every group.
An ordinary fixed-prefix stage keeps its shape when its
literal group and fixed digits change at the same residual budget. The first stage
has a separately checked shape. The verifier tallies
the representative scalar streams of each kind and the
complete two setup/reconstruction pairs, covering the
same register lifetime structure as every other stage.
The largest checked register count is 63.
The literal setup and reconstruction pairs are
[((1004, 3), (1287, 24)), ((1004, 3), (1289, 24))], shown as operation/register-count pairs.
There are 184756 groups per component. Summing both
endpoint pairs and both group streams gives W=174523108758 on every
received word. The checker verifies W^100000 <= 2^3734463,
establishing L=3734463/100000.

## Declared computations

The four runs verify all 4096 tag spectra by four
contiguous ranges [0,960), [960,1960), [1960,3040),
[3040,4096). The arrays store zlib-compressed unsigned
three-byte counts wrapped in Base64. Each run also
verifies its disjoint exact quotient rank sector in a
separate standard-library process, whose successful
exit is checked. The respective support profile sets
are 11–13, 14–16, 17–18 and 19–20. The last run combines
all checked profiles and rank counts to verify the error
bound. It checks the field, group table and all fixed
memory ranges, the operation shapes and exact worst-case
cost exponent. The full read-only banks and the repeated
billion-instruction stream are specified by their
generators and are not materialized by the checker.
