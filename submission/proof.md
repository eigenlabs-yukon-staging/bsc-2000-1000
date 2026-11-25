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

The first 2 tags form the fixed list for that row.
For a received side with quotient q and shift z, add
z to each listed tag d to get the absolute tag t.
Linear J gives s xor J(t)=s0(q) xor J(d), so the list
records give the exact minimum distances for those tags.
The tie rule is fixed and does not change any equal-weight
shell count.

For each ten-element subset of the twenty side slots in
lexicographic order, run through the Cartesian product
of the ten chosen side lists in rank order. Interpolate
the polynomial from those ten selected tags, evaluate
all twenty tags and score the sum of the twenty local
leader weights. Choose the first candidate of minimum
score. Correct the twenty received private words with
the selected side leaders and output the first ten
polynomial values as the tag parameters. The lists are
nonempty, so there is always a candidate. For any fixed
tag tuple this rule finds its nearest codeword in the
private coordinates by independent side minimization.
Every tuple with its tag present in at least ten
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
accounted for, with no later weights. With R=2
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

One component is a binary linear [1000,500] code C.
Consider error patterns whose twenty side patterns are
chosen private leaders and whose ranks fall in the
shortlist on at least ten sides. The set of these
patterns of weight at most t has exact cardinality V_t
given by cumulative coefficients of

    sum_(i=10..20) binom(20,i) A(X)^i U(X)^(20-i).

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
its exact upper rational score 969315802566548673/1000000000000000000.

## Counted 64-bit RAM schedule

The literal tuples emitted by schedule() specify the
decoder. A string denotes a register, an integer denotes
a literal, and move is a free copy. Every other tuple
denotes one named primitive. All branches, loops and
function calls in the generator operate on fixed code
parameters and literal indices. They do not appear in
the emitted instruction stream.

Two read-only field lookup banks begin at 2^24 and
2^25 with 2^24 words each. For each of them combine
the coefficient c and tag t as c*4096+t. The first
bank returns J(c*t), where * is field multiplication,
and the second returns c*t itself. Both bases are aligned
to 2^24. In either bank ORing a tag into its aligned
coefficient subrow gives exactly the stated address.

A read-only group bank begins at 2^26. The ten-element
subsets are in lexicographic order. Each group has ten
21-word rows, one per selected side. The first word
is its side's working shortlist address. The remaining
twenty words are the literal subrow addresses in the
second field bank for the Lagrange coefficients that
take that selected tag to each of the twenty values.
The coefficient generator computes the Lagrange formula,
so every group value is defined. The verifier checks
the table extent, alignment and its address formulas.

The private parity table starts at 2^30 and contains
P(u) with the marker bit for 2^31 set. The subsequent
tag table contains J(t). Two coordinate tables store
the fixed split of the low sixteen and upper fifteen
bits of a side syndrome in the split basis. The quotient
bank has two records per quotient, containing relative
tags below bit 12 and their leader weights above. The
writable arrays store the received private words, marked
syndromes, absolute shortlists, weights and selected-group
contributions. Two words store the incumbent distance
and key across groups. The read-only leader bank at 2^31
has one word for each private syndrome, with its minimum
weight in the low byte and chosen 19-bit private
correction above it. Its exclusive end is 2^32.
The entire fixed layout is disjoint and verified.

For each component load and stitch input word pieces to
form each of its twenty received sides. Separate its
private and parity bits and store the received private
word. XOR the parity-table entry into the observed
parity to get the marked syndrome. Store it. Retrieve
each linear split from the two coordinate tables and
XOR their parts. The high segment contains the marker
in position 15, so the emitted address offset subtracts
2^15. The resulting quotient selects a row and the
tag shift translates its two relative records to the
absolute tags. Write those records to their working list.

For one group, pack its ten complementary marked
syndromes into five 64-bit pairs. For each group row and
each of its two ranks, obtain the absolute tag and
record its weight in reserved working memory. For
each complementary pair, use the group's two literal
coefficients and that tag to read J of the two field
contributions from the first bank. Pack the two 31-bit
displacements into a word with lane stride 32. Store
it for reuse across that group's candidate ranks. This
preparation has four operations per marked syndrome pair.
For each of the twenty selected side/rank records it
has four setup operations and seven operations per
displacement pair. Its total is 800 operations.

A fixed depth-first expansion shares common prefixes of
the ten ranks. Each node loads its next selected weight
and its five packed contribution words. At depth zero,
copy the weight and XOR the five words with the packed
marked syndromes. At each further depth add the weight
and XOR the words into the preceding prefix. For depths
0 through 2, write the six accumulated values into
that depth's reserved scratch slots. Load those six
values when they are needed as the parent prefix.
The other prefix values use named registers. Depth is
fixed at code generation, and the descendants of a
prefix do not modify the saved slots or names of their
ancestors. The leaf then holds all ten marked leader
addresses and the sum of the ten selected weights.

For each of five packed words, split the two addresses,
load both leaders and add both words, six operations per
pair.
A leader word has the exact local weight in its low byte
and the private correction in the higher bits. Add each
unmasked word to the known-side weight sum. The checked
quotient rank spectra cover all local syndromes by weight 11,
so the total weight is at most 20*11=220<256.
All additions preserve the weight sum modulo 256. One
final byte mask therefore yields the exact distance.

At the first leaf initialize the incumbent using free
copies. Every later leaf uses an unsigned strict
comparison and two conditional selects. The first
of equal distances therefore wins. The key is
g*2^10+r, where g is the group and r is the binary
encoding of the ten literal ranks. Each group stores
its incumbent, and every later group loads the two
prior words.

At depth d there are 2^(d+1) nodes. Each uses six
loads for its next record, five XORs, and one weight
addition when d>0. There are another six loads if the
parent prefix was spilled and another six stores if
the node's own prefix is spilled. In addition each
leaf has 5*6+1+3 operations, except the first which
skips the three incumbent comparisons and updates.
The first group therefore costs 60417 with its
setup and last two stores; later groups cost
60422 with setup, last stores and the two
loads of the prior best values.

After the last group shift the winning key by ten
to recover its group and extract the individual ranks
from the low bits. From its group bank rows read the
chosen absolute tag for each selected side, through
the stored working list pointer. Read each of the
twenty Lagrange contributions from the second field
bank and XOR together the ten contributions at each
side. Those are the same twenty tags that were scored.
For each side look up J of its tag, XOR it with the
marked syndrome, load its leader, extract the private
correction and XOR it into the received private word.
This is the component decoder of the error proof.

Pack the twenty private words and first ten tags into
successive fields of 19 and 12 bits, respectively.
A field start on a word boundary uses a free copy.
A field crossing a word boundary uses shifts and OR
to finish and store the preceding word and keep the
high part for the next word. The first component stores
its partial output word and the second reloads it at
the shared boundary. Its own last partial output word
has only 40 used low bits, so all 24 padding bits are zero.

All stages use the same fixed generation-time register
allocator. Treat each new destination as a fresh value,
record its last operand use, read operands before assigning
the result and return registers after their last use.
This retains the correct virtual value for every operand.
Use separate allocations for component preparation, each
group, and component reconstruction. Their carried data
is explicitly written to and read from working memory.
The register lifetime shape is identical in every group
after the first, since changing the group changes only
integer literals, never operation names, virtual operands
or the count and order of the rank loops. The first group
has its own fixed shape. The checked register totals are
50 in group scoring and no more than 24 in the two
component endpoint stages. All operands fit 64-bit words.

The computation checks those group formulas against
literal representative schedules of both shapes; the
liveness argument makes this exact for all groups.
It generates the full setup and reconstruction streams
for both components and counts their operations separately.
Their counts are (604,1287) and (604,1289), respectively.
Thus the fixed decoder executes, on every received word,

    W=604+1287+604+1289+2*(60417+184755*60422)=22326657838.

The checked comparison W^20000 <= 2^687561 gives
L=687561/20000.

## Declared computations

The four runs verify all tag spectrum rows, in the
intervals [0,960), [960,1960), [1960,3040), [3040,4096).
The four files store compressed exact coefficient arrays
with unsigned three-byte little-endian coefficients.
Each run also checks one of the four disjoint quotient
rank sectors. That rank check uses a standard-library
child process, in parallel with tag spectra and onward
arithmetic. The parent checks its exit status. Every
run recomputes and verifies its share of the support
profiles in profiles.json. The respective support
size groups are 11–14, 15–17, 18–19, and 20.
The last run combines those profiles and exact rank
rows, checks the occupancy and error bound with integer
arithmetic, verifies the field and memory ranges, and
checks the counted instruction shapes and exponent.
It does not materialize the full table banks or expand
the billions of repeated ordinary group instructions.
