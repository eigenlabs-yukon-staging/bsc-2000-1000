# Twenty-sided polynomial tags with one extended list

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

The first four relative tags are stored for each quotient.
The absolute tags for a received side are obtained by XOR
with the received tag coordinate z. For a listed relative
tag d, s xor J(d xor z)=s0(q) xor J(d), so the record
weight is the exact minimum error weight for that tag.
Select the side with the smallest third entry weight.
Break ties by the larger second entry weight, then larger
fourth entry weight, then smaller side index. Use all three entries there and
only the first two at every other side.

For each ten-element subset of the twenty slots, in
lexicographic order, enumerate the Cartesian product of
its selected side lists in rank order. Interpolate the
degree-at-most-nine polynomial from these ten tags and
score the twenty sides by their exact local leader weights.
Keep the first minimum-score tuple. Apply the chosen
private corrections and output the first ten tag values
as the tag parameters. Every tuple whose tag occurs in its
side lists on at least ten sides is included. For a fixed
tag tuple, minimizing the private distances by the chosen
leaders gives a closest codeword having those tags.

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
counts by weight and rank through weight eleven. The fixed local ranking determines the eligibility tests.

Each quotient sector has exactly 2^29 syndromes and its
twelve checked shell totals sum to 2^29, so all local
minimum-leader weights are covered. Let T(X) record
the number T_w of these leaders at each local weight w,
summing over all four quotient sectors.

## Exact admitted volume and two-error overlap

One component is a binary linear [1000,500] code C.
Consider patterns whose twenty sides are chosen private
leaders and whose ranks pass the actual list eligibility
test on at least ten sides.

### Exact volume with the selected extra rank

For each quotient q let w_i be the weight of its ith
relative tag leader in one-based rank. It is the first
weight with quotient population at least i. The priority
key p(q) is (w3<<10) OR ((31-w2)<<5) OR (31-w4).
Since every weight lies between 0 and 31, the fields are
disjoint and the order of keys is the stated lexicographic
order. The exact rank checker classifies each completed
row by this key. In class p and at weight w its counts
of ranks less than two, rank exactly two, and all ranks
are the respective population increments
min(b,2)-min(a,2), 1[a<3<=b], and b-a for successive
populations a,b. Each sector checks these values directly
against priority.json.
The totals across
classes also match the unclassified rank counts.

An error pattern on a side determines the observed
quotient: adding the true tag displacement leaves the
quotient unchanged. The absolute tag shift also leaves
the weights of the relative list entries unchanged.
Thus the priority used by the decoder is exactly this
error quotient's priority. On a side whose error is a
chosen leader, the true tag is in the prescribed list
if its relative rank is below two, or if it has rank two
and that side receives the extra entry.

For each priority class p define two bivariate
polynomials A_p(X,Z), B_p(X,Z). Both count all its
side leaders by noise weight in X. In A_p a leader
has an extra factor Z exactly when its rank is below
two. In B_p it has Z exactly when its rank is below
three. For any fixed positions occupied by class p,
their product is A_p^n when none receives the extra
entry, or B_p*A_p^(n-1) when one receives it. There
is no additional factor n in the second case, because
the selected position is the first of those positions
by side order. All positions have the same shell
counts, so any particular choice of that position
would give the same product.

For a selected side position i with priority class p,
every earlier position must have priority above p and
every later position must have priority at least p.
Let U_p be the sum of plain polynomials A_q over q>p.
Then summing the count over the twenty possible positions
gives the class polynomial

    B_p * sum_(i=0..19) U_p^i (U_p+A_p)^(19-i).

Every complete assignment picks exactly one selected
position and class, so summing these polynomials over p
counts all assignments and all local chosen-leader
choices once, with the true number of admitted sides.
Use the cumulative coefficients through X-degree t with
Z-degree at least ten for V_t.

Multiply the finite geometric sum by A_p to obtain
(U_p+A_p)^20 - U_p^20. The verifier traverses classes
from largest to smallest, maintains U_p, and evaluates
that identity in its packed integer representation.
It checks the integer quotient is exact by multiplication
before multiplying by B_p. The geometric sum itself
has nonnegative coefficients and is within the same
radix bound as the other counts. This rescales no
coefficient and uses no approximation.

The exact checker stores the exponent X^w*Z^k at
index w+221*k because the total weight is at most
20*11=220. Polynomial arithmetic uses base 10^190,
with coefficients below it since even all 20 choices
of 2^31 local leaders give fewer than 10^190 terms.
The implementation uses standard-library decimal
integers with sufficient precision for all product
digits, trapping any rounding or inexact result.
It reads the final radix blocks as exact integer
coefficients and checks their total equals (2^31)^20.
Thus the calculation uses no floating-point estimates.

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
difference is identically zero, all twenty patterns are
equal, which is excluded. The tag difference otherwise
comes from a nonzero degree-at-most-nine polynomial;
by the root bound it is nonzero on at least eleven sides.
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
its exact upper rational score 965525619774802117/1000000000000000000.

## Counted scalar RAM schedule

The emitted instructions use one tuple per named RAM
primitive, except that move is a free copy and labels
mark program locations. Literals, including resolved
branch targets, are fixed by the generator. The code
uses divu and modu for unsigned division and remainder.
All loops in the generator are expanded in the emitted
program. The only data-dependent control is the
explicit conditional branch dispatch below.

The two field banks beginning at 2^24 and 2^25 each
contain 2^24 words. For a coefficient c and tag t the
first stores J(c*t) and the second stores c*t at
offset c*4096+t. The upper bank begins above the
lower bank, and ORing t into a coefficient row's
aligned start forms exactly the intended address.

The group bank begins at 2^26. Each ten-side group
has ten 21-word rows. A row's first word gives the
working list address for that side. Its remaining
twenty words give the second field bank's aligned
coefficient row addresses for the Lagrange polynomial
of that side at the 20 points. The field formulas
determine every word. The remainder of the bank
layout is the marked parity table at 2^30, tag and
linear split tables, and quotient lists of four
records per quotient. A record has a relative tag
in its low twelve bits and its exact weight above.
The writable arrays give separate addresses for
private data, marked syndromes, absolute lists,
group weights and packed contributions, spilled
prefix values, incumbent distance and key, and
the selected extra side index. The leader bank at
2^31 contains 2^31 words, with exclusive end 2^32.
A leader word stores its exact weight in the low
byte and the nineteen-bit private correction above.
The disjoint address bounds are checked exactly.

At a component start, assemble each 50-bit received
side and store its private portion. XOR P of that
portion and the 2^31 marker into its received parity
to obtain and store the marked syndrome. Split it
by the two linear table lookups; the high split
address subtracts 2^15 to compensate for the marker.
Read the four relative list records and XOR the
tag shift into their tag fields. Keep all four
absolute records for the priority decision.

Read the second, third, and fourth weights at each side
and form the packed priority key above by scalar shifts,
subtractions and OR operations. Each side costs twelve
operations for these reads and key formation. The first
key and side initialize the minimum by copies. Each
later comparison uses an unsigned strict less-than and
two selects, so equal keys keep the smaller side index.
Store the chosen index. This selection is 298 operations.
After a load of the index, test it for equality with
0 through 18 in order and branch to the corresponding
case on a match, or fall through to case 19. Only that
case's group stream is executed. Each case finishes
by jumping to the shared reconstruction label. Thus
case i consumes 2+2*min(i+1,19) dispatch operations,
counting the load and final jump. Each case has the
same fixed set of groups and fixed ranks determined
by its literal extra side, without runtime calls.

For each group pack the ten complementary marked
syndromes into five word pairs with half stride 32.
For every eligible selected-side/rank entry, read
and store its weight and form five packed displacement
contributions through the first field bank. Each pair
needs two address operations, two loads, a shift, an OR
and a store. If position d has r_d eligible ranks, the
setup cost is
20+39*sum_d r_d, including four operations per marked
syndrome pair. At most one r_d is three; all others
are two. If the extra side is absent all are two.

The case stream expands the Cartesian product of
its ranks in lexicographic order, sharing contributions
at each prefix. At depth zero each node loads its
next weight and five words, copies its weight and
XORs the five words into the packed marked syndromes.
At each later depth its six loads obtain the next
entries, its five XORs add the displacements to the
prefix, and its addition extends the weight sum.
Depth zero stores its six prefix words to fixed slots.
A child of depth zero loads those six words to recover
its parent prefix.
Descendants cannot overwrite their ancestors' slots.
All other prefixes have distinct virtual registers
by depth. The number of nodes at depth d is
the product of r_i for 0<=i<=d.

At a leaf split the five words into the ten marked
leader addresses. Load and add the two leader words per
pair, for six operations per pair including the split.
A leader word has the exact local weight in its low byte
and the private correction in the higher bits. Add each
unmasked word to the known-side weight sum. The checked
quotient rank spectra cover all local syndromes by weight 11,
so the total weight is at most 20*11=220<256.
All additions preserve the weight sum modulo 256. One
final byte mask therefore yields the exact distance.
The first leaf assigns the incumbent via
free copies. Every other leaf uses a strict comparison
and two selects, so a tie preserves the first candidate.
The winner key is g*3^10+o, where g is its group and
o is the fixed ten-digit base-three encoding of its
rank tuple. This encoding is injective even for groups
with fewer than three ranks in some positions.
Each group stores the two incumbent words and all
but the first group load them before the search.
The first group is five operations cheaper than an
ordinary group with the same radices (two loads and
three first-leaf update operations).

For a case with extra actual index i, there are
binom(19,10) groups without it. The number with
the extra at position d is
binom(i,d)*binom(19-i,9-d), with out-of-range terms
zero. The verifier checks the node, setup and leaf
counts and generates one representative of each
first and ordinary group shape. The operands' lifetime
shape depends only on those radices and on whether
the group is first, so the register bound applies
to every actual group. Sum the ordinary costs with
these multiplicities, subtract five, add the exact
dispatch cost and take the maximum over i. The
component repeats its independently bounded preparation
and reconstruction before and after these operations.

Read the winning group and ranks by unsigned division
and remainder, fetch its selected absolute tags and
apply the second field bank contributions to recover
all twenty tags. Use all winning tags and their marked
syndromes to read the chosen private corrections.
Output twenty corrected private words and the first ten
tag values in nineteen- and twelve-bit fields.
Store any partial output word at the component boundary
and explicitly load it in the next component. The final
word has 40 data bits and 24 upper zeros.

The scalar stream within each setup, selection, group
and reconstruction receives a separate generation-time
physical register allocation. Number virtual destinations,
record their last use, and release an operand's register
after its last read. Reading before each result write
makes reuse safe, so induction gives exactly the
virtual program's values. The carry between stages uses
only the explicitly stored working data. Dispatch itself
uses two registers. Representative groups of all shapes
and all complete endpoint/selection stages are checked,
giving a maximum of 62 registers.

The checked worst-case group-and-dispatch maximum for
one component is 13870741455. Setup and reconstruction
counts are (764, 1287) and (764, 1289) and the selection is 298
per component. Summing these endpoint and decision
costs gives worst-case W at most 27741487610.
The exact inequality W^100000 <= 2^3469133 establishes
L=3469133/100000.

## Declared computations

The four runs check the full 4096 local tag spectra
in intervals [0,1100), [1100,2200), [2200,3350),
[3350,4096). The spectrum arrays are unsigned
three-byte coefficients compressed with zlib and
wrapped in Base64. In parallel each run recomputes
one disjoint quotient sector, with the total rank
histogram and its priority-conditioned shell counts.
The parent checks the child's successful exit.
They also verify support profiles 11–13, 14–16,
17–18 and 19–20, respectively. The last run combines
the profiles and class counts to check the exact
admitted volume, error and worst-case operation
claims, including field and memory formulas and
the fixed operation and register shapes. It does
not materialize the full preloaded banks.
