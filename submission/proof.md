# Two cyclic fields with exact weight-indexed admission

## Encoder and decoder

Use the polynomial basis of F2[X]/(P), where
P=X^500+X^21+X^13+X^8+1. The declared field check applies
the Rabin criterion with prime divisors 2 and 5 of 500.
Split the message into two 500-bit field words b0 and b1.
Transmit both, then the parity fields

    p0=L(b0)+b1^3,   p1=L(b1)+b0^3.

The encoder is injective by its 1000-bit systematic prefix,
and use of all four 500-bit fields gives exactly 2000 bits.

In each field, L is block diagonal with twenty 24-bit cells,
one 12-bit cell and two 4-bit cells in increasing bit order.
For a 24-bit input u use the binary polynomial G=9205487
(in low-coefficient-first order), of degree 23. The low 23
bits of its parity are r=(u*X^23) mod G. The high parity bit
is the parity of wt(u)+wt(r). For the 12-bit Golay cell use
the same rule with degree 11 and polynomial
1+X^2+X^4+X^5+X^6+X^10+X^11. For each 4-bit Hamming cell
the parity columns are the four words 14,13,11,7. All are
linear. Construction.py gives the precise integer ordering.

The decoder enumerates all anchors x0 of distance d<=100
from received Y0. For each use entry C(d) of construction.py BUDGETS as the initial remaining-cost state. Starting from the first child cell,
form its systematic and adjusted parity observations u,v
using received Y1 and P1+x0^3. At width m define the local
syndrome s=v+L_m(u). For each error e on the systematic
part, the paired distance and unnormalized weight are

    d_s(e)=wt(e)+wt(L_m(e)+s),
    A_s(e)=9^(2m-d_s(e)),   S_s=sum_e A_s(e).

Assign cost c_s(e)=ceil(16 log2(S_s/A_s(e))). The preloaded
integer cost is bit_length((S_s^16-1)//A_s(e)^16), since
the quotient is below 2^c exactly when S_s^16<=2^c*A_s(e)^16.
Costs are positive because every complete row has more than
one positive term.

All candidates of widths 4 and 12 appear in fixed rows
ordered by (distance, numerical systematic error). In width 24
use the complete row of 2^24 noise patterns, through an explicit
syndrome orbit map and coordinate permutation defined below.
This preserves the distances and the syndrome. The selected
cell is u+e, where e is the systematic part of its noise pattern
in original coordinates. All rows are in nondecreasing distance
and thus nondecreasing cost.

Depth-first search the 23 child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(23-i), where r is
the current state and p is the syndrome profile type for this cell.
The child's state is E_(i,p)(r) minus the selected cost.
Each retained entry has positive integer cost. The reserve test
allows the implementation to stop at the first failing entry in cost order.
Score each complete traversal by full received distance. Start with the
all-zero codeword as fallback and replace the incumbent only on a strict
improvement in the stated traversal order. The proof below derives the
table entries for every possible state.

## Extended 48-bit cell data and its orbits

For this cell use an auxiliary raw 48-bit coordinate order:
low 23 parity bits, 24 systematic bits, then the overall parity
bit at position 47. The codeword body is divisible by G and
its full parity is even. The 24 shifted copies G*X^i for i=0..23,
extended by their overall parity bit, lie in the code. They
are independent by their distinct leading degrees, hence a
basis because the code has dimension 24 from its systematic
definition. The declared computation checks all pairwise binary
dot products are zero, so the linear code equals its dual.

The coordinate permutation group consists of matrices over F47
of determinant one, modulo multiplication by -1, acting by
x -> (a*x+b)/(c*x+d) on F47 and the point at infinity 47.
A nonzero denominator uses its field inverse, zero denominator
maps to infinity, and the image of infinity is a/c for c nonzero
and infinity otherwise. The group list is the deterministic
loop order in qr.py, keeping the lexicographically earlier of
each opposite pair. It has 51888 elements, enumerated exactly.

Translation x->x+1 and inversion x->-1/x preserve the 24
basis codewords' span, as explicitly checked. These two
permutations generate the group: their products give every
upper triangular elementary matrix U(t) with diagonal ones,
and conjugation by the inversion matrix S gives the lower
elementary matrices V(t). These eliminate the off-diagonal
entries of any determinant-one matrix, after applying S if its
upper-left entry is zero. The remaining diagonal matrix with
entries a and a^-1 equals U(a)*V(-1/a)*U(a)*S. Thus every
listed permutation preserves the code.

For the syndrome s in the actual decoder, each low parity
coordinate and the overall parity coordinate give a different
unit syndrome, since the systematic bits are zero there.
The syndrome of any raw error is the XOR of its columns
specified by parity(e_sys)+e_par. On these 24 pivot error bits,
the syndrome action of any preserving permutation is therefore
the binary map whose columns are the syndromes of the permuted
pivots. It is well-defined on every coset because the code is
preserved. Its set of fixed syndromes has size 2^(24-rank(M-I)).
The checker enumerates all matrices, calculates these ranks
by binary elimination and uses Burnside's orbit count. The
result is 426. Burnside's formula follows by counting pairs
(s,g) with g(s)=s: summing over s in any orbit contributes
the group size.

For each syndrome orbit its representative is the smallest
numerical syndrome and orbit indices rank representatives in
increasing order. For a syndrome s with representative r, store
the orbit index and the first group element sending r to s.
The data generator scans syndromes in increasing order, labels
the image of each as-yet-unlabelled representative under all
permutations in group order, and makes every assignment once.
Thus this fixes the mapping for every s and gives all 426 rows.

For each representative, the fixed row lists the 2^24 raw
noise patterns with that syndrome in increasing (distance,
numeric raw noise) order. For each systematic error e, its
parity error is L_m(e)+r, so the full list and ordering are
explicit. Store the systematic error in 24 bits followed
by its six-bit total distance. Two consecutive 30-bit
entries are packed into each word, with the even rank in
the low part. The representative and the parity table
recover the complete raw pattern. Every orbit has a
complete weight profile, and entry cost at distance d is
computed from that profile's normalizer and 9^(48-d).
The separate cost table stores these integers. The runtime
maps the recovered raw noise pattern back to syndrome s
using the stored coordinate permutation. That permutation
preserves weight and is a bijection, so every coset weight
profile is realized without truncation.

## Exact local syndrome spectra

Let a range through the 24-bit labels in the displayed code
basis, and let v(a) be the corresponding 48-bit codeword.
If a label is used as a generator of dual constraints, the
linear syndrome of a raw error e has entries its dot products
with the basis rows. This syndrome differs from the actual
systematic syndrome by an invertible change of basis since
both have kernel the same code. Hence their multisets of
coset weight enumerators are identical.

For a dual syndrome label s and each weight d, the number of
weight-d errors in its coset is

    h_s(d) = 2^(-24) * sum_a (-1)^(a·s) K_d(wt(v(a))),
    K_d(w)=sum_j (-1)^j binom(w,j)*binom(48-w,d-j).

Indeed binary character orthogonality tests each of the
24 parity constraints. The signed sum over weight-d errors
for a fixed dual word counts the intersection weight j
with its support, which gives K_d(w). The checker enumerates
all 2^24 codeword weights as XORs of the independent lists for
two 12-generator halves. It finds only weights
0,12,16,20,24,28,32,36,48.

The signed transforms of the weight indicators for 12,16,20
and 24 are enough. The weight-zero label is zero; the weight-48
label is the unique quotient of the 47 all-one body bits by G,
extended to even parity. Adding that codeword complements
coordinates, so the transform for weight 48-w equals the one
for w times the sign of this label against s. The transform
for weight 48 alone is exactly that sign.

For an indicator row, the experiment computes the complete
Walsh--Hadamard transform by successive butterflies, with
one radix-2^32 digit per input label. It stores each signed
integer as its residue modulo 2^31. At each stage collect the
low and shifted high half of every active butterfly by a mask.
The low output is their sum, reduced modulo 2^31 by the
same mask. The high output is their difference; before integer
subtraction put a guard bit 2^31 in every active low digit to
ensure no lane borrows. Mask again and shift back to the high
positions. A sum is below 2^32 in each digit. Inactive digit
positions are zero. This is exact packing of all scalar
butterflies. Every transform prefix has magnitude at most 2^24,
so residues determine signs and values without ambiguity.

Group the exact tuples of four computed transform values
and the all-one-label sign. For each tuple generate all 49
weights of the coset enumerator by the displayed formula,
checking nonnegativity and exact division by 2^24. The first
experiment reproduces each of the 14 profile types with its
exact multiplicity in profiles.json. Later experiments consume
that same frozen data. No estimator or sampled row is used.

## Child admission by exact noise weight

On the true path the predecessor cube equals the true
cube. Subtracting it leaves the original BSC(1/10)
observations. Every cell syndrome on that path depends on
its own noise, independently of previous true cell noise
and the anchor. Cell i has b_i observed bits, with b_i=2m
for a paired width-m cell and b_i=1 for a systematic-only
cell. Let H_(i,p,c,w) be the number of retained true noise
patterns of weight w and nominal cost c across all
syndromes of its profile type p. Multiply the retained
count for one syndrome by its profile multiplicity.
The QR orbit permutation preserves total weight and cost.
In a systematic-only cell there is one type, and one
pattern each at weights zero and one with their stated
costs. The channel mass for a term (p,c,w) is

    m_(i,p,c,w) = H_(i,p,c,w)*9^(b_i-w)/10^b_i.

For each i, use the same effective cap E_(i,p)(r) and
remaining state r as the decoder. With s=E_(i,p)(r),
a true pattern of cost c passes precisely when it is retained in
its row and c<=s-(N-i), where N is the fixed number of
child cells in the program. The next state is s-c.
On the true path these tests apply directly,
since the predecessor adjustment and the local noise
syndrome coincide with their actual values. Define
P_i(r;X) with coefficient at X^z equal to the probability
that the suffix from i is admitted and has true noise
weight z. The exact recurrence is

    P_(N+1)(r;X)=1
    P_i(r;X)=sum_(p,c,w: c<=E_(i,p)(r)-(N-i))
                   m_(i,p,c,w) * X^w * P_(i+1)(E_(i,p)(r)-c;X).

This conditional law depends solely on r and the remaining
independent true noise. The prefix reserve is part of every
transition, and unlisted true patterns have no term.

The computation derives lower integer coefficients at
scale S=2^64. Round each local mass down to
floor(S*m_(i,p,c,w)). The empty suffix polynomial has
coefficient S at degree zero. At a step multiply a lower
suffix coefficient by the appropriate local mass integer
and sum contributions, then floor the resulting coefficient
divided by S. Nonnegativity gives, by induction, lower
bounds on S times all the exact P coefficients. Only degrees
through 230 are required, because the scoring contribution
in this proof only uses total noise weights up to 230.

The verifier stores each polynomial with base 2^136, using
one digit per weight. Multiplication by the local mass
integer and a shift by 136*w combines the two factors in
the suffix recurrence. Before dividing at a stage, each
coefficient is at most S^2, since the total local mass is
at most one and every lower suffix coefficient is at most
S. All contributions are nonnegative. There are therefore
no carries among 136-bit digits. Discard digits above the
cutoff, shift the whole integer right by 64, and mask the
low 65 bits of each digit at its original position. This
equals flooring each original digit divided by S: displaced
bits from a following digit begin at bit 136-64=72 and
are erased by the mask. This implements the rounded
probability recurrence using only exact integer operations.

For a true anchor depth d within the radius, start the suffix
at C(d). Let l_(d,z) be its lower integer probability numerator
at child weight z, with denominator S=2^64.
All child noise patterns of weight z on n_c observed bits
have the same probability 9^(n_c-z)/10^n_c. Thus at least

    ceil(l_(d,z)*10^n_c/(S*9^(n_c-z)))

patterns of that weight pass child admission. Rounding up is
valid because the actual count is an integer above the scaled
lower bound. Subtract this count from binom(n_c,z) for an
upper bound on the number that fail. For depth d there are
binom(ANCHOR_BITS,d) anchor patterns, and the unaffected anchor
parity contributes binom(M,v-d-z) patterns at total weight v.
The verifier sums the product over d within the radius and
z<=v-d. This gives an integer upper count of removed true
noise patterns at each total weight v through 230.
Part 2 computes and checks all the resulting counts in misses.json.

## Competitors changing both fields

Write q=2^500 and Q=q^2. For nonzero e, the cube
difference on a uniform x is

    (x+e)^3+x^3 = e^3*((x/e)^2+x/e+1).

The F2-linear map z->z^2+z has kernel {0,1} and
image the kernel of absolute trace. Trace is
F2-valued and nonzero by the polynomial root bound,
and its value on z^2+z is zero. Hence the difference
is uniform on one affine hyperplane with mass 2/q
on each of its members. Its normal in the binary
basis has entries Tr(e^-3 X^j). The map from a
field multiplier to the normal is injective since
a nonzero multiplier permutes the field and cannot
annihilate trace on a spanning set. For each
nonzero normal mask there are at most three
directions e giving it by the cubic root bound.

A full-change competitor is specified by nonzero
directions e_i in both fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the two translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 1000
parity bits can be evaluated by expanding the
two affine indicators. Multiplication by their
common point mass 4/Q leaves 1/Q times a sum
of four terms. The empty mask term is the
parity ball volume. Each nonempty term has an
absolute bound equal to the absolute signed sum
of the parity ball against the combined mask of
its selected parity components.

For a 1000-bit mask of weight w>0 the signed ball
sum is, up to translation sign, the coefficient
K_l(w) of X^l in (1-X)^(w-1)*(1+X)^(1000-w).
Indeed the polynomial of shell signs has the same
form with exponent w on (1-X), and summing the
prefix divides by 1-X. Let N=999 and
g=N-2*(w-1). Differentiating this polynomial F gives
(1-X^2)F'=(g-NX)F. The coefficient recurrence is

    (l+1)*K_(l+1)=g*K_l-(N-l+1)*K_(l-1),

with K_0=1 and K_(-1)=0. The calculation checks each
division exactly.

In a term that selects k masks, their combined
weight w has at most the following multiplicity
across ordered directions on the contributing
systematic fields,

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(500*h,w).

The alternating sum counts weight-w combinations
with no zero mask on any of the k blocks. This is
inclusion and exclusion from all binary weight-w
words on the k blocks. Each mask has at most three
direction preimages. For a given l and k define
H_k(l,n) to be the largest sum of n absolute
magnitudes in this list, sorted by |K_l(w)| with
the displayed multiplicities. The computation takes
full groups of largest magnitude then a partial
group, giving the maximum over any at most n
directions. The available multiplicities exceed the
capacities used.

Fix the true noise. In a chosen k-mask term the
relevant k systematic fields are the predecessors
of the selected parity fields. For a combined
systematic error weight a on them, there are at
most binom(500*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(500*(2-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(500*(2-k),b) *
    H_k(t-a-b, binom(500*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(2,k) selections of masks, and add
the 2000-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_fixed(t). The expected
number of full-change competitors within total
t for each fixed true noise is at most D_fixed(t)/Q.
The checker verifies D_fixed(t) is nondecreasing through t=230.

## Averaged full-change ball intersections

The occupancy calculation needs a bound on the total count
of ordered full-change pairs having a common received word
within distance t of both codewords. At cutoffs 212 through
230 we also bound this count by summing intersection sizes,
which averages over all true noise patterns in the ball.

Write M=500 and K=1000 for the lengths of one field and
the complete systematic block. Let V_t be the 2000-bit
Hamming ball volume. For two binary words at distance d,
the intersection of their radius-t balls has size J_t(d).
On their d differing positions, let a of the received bits
differ from the first word. The other d-a differ from the
second word. The unaffected positions can contribute at
most t-max(a,d-a) differences. Thus

    J_t(d) = sum_(a=0..d) binom(d,a) *
                  sum_(b=0..t-max(a,d-a)) binom(2000-d,b).

Use an empty inner sum for a negative bound, so J_t(d)=0
for d>2t. Every J_t(d) lies between zero and V_t.

Fix both nonzero systematic changes. Write d_s for their
combined weight. Their parity difference law, averaged over
the true message, has the uniform term 2^-K and one signed
character term of that magnitude for each nonempty choice
of supported mask blocks. This is the same derivative
expansion as above. Any term with a selected mask of total
weight w contributes, in absolute value before dividing
by 2^K, at most

    | sum_(p=0..2t-d_s) J_t(d_s+p) K_p(w) |,

where K_p(w) is the coefficient of X^p in
(1-X)^w*(1+X)^(K-w). This follows by summing the character of
weight w over each parity shell of weight p.

For j=min(w,K-w), pair j factors (1-X)(1+X).
Taking the absolute value in each coefficient bounds the
coefficients of the signed polynomial by those of

    U_j(X) = (1+X^2)^j * (1+X)^(K-2j).

Each increase of j replaces one factor (1+X)^2 by
(1+X^2), so U_j decreases coefficientwise. Round j down
to its bin endpoint in the fixed sorted MASK_BINS. The
nonnegative convolution of J_t(d_s+p) with the appropriate
U_j coefficients therefore bounds the term for every
mask weight in that bin. The common upper value is
nonincreasing as the bin endpoint increases.

Use the same rank argument as for the per-noise bound.
For one nonzero mask block at weight w, its multiplicity
among first directions is at most 3*binom(M,w). For two
selected blocks at combined weight w the capacity is at
most 9*(binom(K,w)-2*binom(M,w)), since both masks are
nonzero. Sum these capacities over each bin. When only one
block is selected, for the unconstrained systematic change
of weight b there are binom(M,b) choices. Among the other
systematic changes there are at most binom(M,d_s-b)
choices. Use that many top signed-term bounds from the
one-block capacity list, then multiply by binom(M,b) and
sum over b with both weights positive. Multiply by two
for selection of either parity block. When both blocks
are selected use the top binom(K,d_s)-2*binom(M,d_s)
bounds from the two-block list, the number of systematic
changes with both blocks nonzero. Sum these corrections
over d_s=2 through 2t to obtain R_t. These counts ignore
signs and take the largest allowed terms, an upper bound
for any phases arising from L and the derivative.

The packed convolution uses reversed J_t coefficients,
so its coefficient at degree 2t-d_s equals the required
nonnegative sum in U_j. Both constituent polynomials for
U_j have nonnegative coefficients, whose product sums to
2^(K-j)<=2^K. Each coefficient after multiplying by
reversed J_t is at most V_t*2^K. The fixed digit stride
exceeds its bit length, preventing carries in either
multiplication. Truncating above degree 2t preserves
exactly the coefficients used. The ranking uses integer
cumulative counts and exact top partial sums.

The empty mask contribution, summed over all systematic
differences including the omitted zero-block directions,
is at most V_t^2. To see this, sum J_t(wt(z)) over all
2000-bit differences z. Each ordered pair of vectors
in the ball contributes once at its XOR difference. The
parity density factor 2^-K cancels the count Q=2^K of
true messages. Consequently the total ordered full-change
pair count is at most V_t^2+R_t. It is also at most
V_t*D_fixed(t) by the per-noise calculation.

At the specified cutoffs take the smaller of D_fixed(t)
and ceil((V_t^2+R_t)/V_t). At other cutoffs use
D_fixed(t). Finally take the running maximum of these
values through t and call it D_t, so it is nondecreasing.
For every t the ordered full-change pair count is bounded
by V_t*D_t. The experiments derive all nineteen integer
ceilings in pairs.json. Parts 1 and 2 check cutoffs
212 through 224 and 225 through 230, respectively.

## Competitors changing a proper subset

Let q=2^500 and Q=q^2. For a proper nonempty set of changed
systematic fields, its cyclic runs of consecutive changes have disjoint
affected coordinates. A run (i,\ldots,i+a-1) affects its (a)
systematic fields and the (a+1) parity fields (i,\ldots,i+a).
The distance change of the whole competitor is the sum of the changes
on these disjoint sets. Thus any proper competitor of distance at
most that of a given true word implies a single-run competitor of
distance at most that of that word. We bound removal for inclusive
comparisons, including all ties.

Fix the true noise and a run of length a=1.
In each cell of the first affected systematic and parity pair, their
joint syndrome (s=e^{\rm par}+L_m(e^{\rm sys})) is invariant under
this run: their differences are (p) and (L_m(p)) for the candidate
direction (p), since its predecessor is unchanged. The possible
weights for the resulting error on that pair have polynomial
(H_s(Z)=\sum_d h_s(d) Z^d). For each of the other (a) affected parity
fields, the fixed-direction parity difference is an affine uniform
hyperplane variable with mass (2/q) at each supported output.
These variables are independent over the uniform transmitted word
because they use distinct true predecessor fields. The remaining
(a-1) affected systematic fields and those (a) parity fields
jointly have (u=500(2a-1)) bits. Each potential systematic error
pattern fixes the direction on the run. We can upper bound the
expected number of competitors of local distance at most (w) by
((2/q)^a) times the coefficient prefix through (w) of

    (1+Z)^u * product_(23 pair cells) H_s(Z).

Here (w) is the local true noise weight in all (1000+u) affected
bits. Allowing zero directions in the polynomial only adds
nonnegative terms. This use of the full set of error strings never
assigns a hyperplane point mass to an actual direction equal to zero:
the bound (2^a/q^a) applies to the actual run directions (all
nonzero) before enlarging the final counting sum.

For (0<r<1), evaluate the polynomial at (Z=r^2). Each term of
degree at most (w) has a factor at least (r^{2w}), so the prefix
is at most (r^{-2w}(1+r^2)^u\prod H_s(r^2)). The probability
of existence is at most the lesser of 1 and the resulting expected
count. The inequality `min(1,x) <= sqrt(x)` for (x\ge0)
therefore gives the bound

    2^(-floor(499*a/2)) * r^(-w) * (1+r^2)^(u/2)
        * product_(23 pair cells) sqrt(H_s(r^2)).

This rounds the power of two upward. Use the rational r=42/100.
The verifier bounds each square root from above at resolution
(2^{-64}). Write an evaluated (H_s(r^2)=n/d) and take (j)
to be the integer square root of the floor of (2^{128}n/d).
Then the exact numerator of the rounded bound is (j) if
(j^2d=2^{128}n), or (j+1) otherwise.

For each length, to average over the channel noises, form the
polynomial of each paired cell with coefficient at degree (d)
equal to the sum of (h_s(d)) times the rounded integer root
numerator over all syndromes (s). Multiply the twenty width-24, one width-12 and two
width-4 polynomials and ((1+Z)^u). At local degree (w)
apply (r^{-w}), via a common denominator using the largest local
degree that is needed, 230. Multiply by ((1+Z)^{2000-1000-u})
to include outside weights. A final coefficient at total degree
(v) is multiplied by the common run factors and the channel
factor (9^{2000-v}/10^{2000}). Denote this bound on the joint
run-removal probability at weight (v) by (J_{a,v}). The
translation-invariant channel makes the bound identical at each of
the two starts. Terms above total degree 230 are never needed.

In this computation the rounded root numerator is at most (2^{76}):
(H_s(r^2)\le H_s(1)=2^m) and (m\le24).
The sum of the integer coefficients of the paired factors
is at most (2^{1000+23\cdot76}), and the binomial factor contributes
(2^u\le2^{500}). Thus the first product is safe in 5000-bit digits.
Writing (r=a_0/b_0), the local weight multiplier on the integer
coefficients is (b_0^w a_0^{230-w}\le100^{230}<2^{1529}).
After the remaining binomial factor the resulting sum is bounded by
(2^{2000+23\cdot76+1529}<2^{7000}). Truncations and binary powers
therefore do not create carries in the exact product.

## Occupancy after removals by exact noise weight

At each total noise weight v, upper-bound the mean number
of removed patterns per true message. Remove a pair of true
and received words if the anchor or child admission test fails,
or if some single-run competitor has distance no greater than the true word. Two retained words
that differ on a proper subset cannot coexist, since at least the
farther is removed and both are removed on a tie.

For the anchor, sum

    binom(ANCHOR_BITS,d)*binom(2000-ANCHOR_BITS,v-d)

over d above the radius and at most v. The child count was
computed for the disjoint range within the radius. For the
proper-run term, use the convolution over unaffected noise
coordinates described in the previous section. Before
multiplying by the common channel probability
9^(2000-v)/10^2000 its value is the upper mean count of
removed patterns at weight v. Sum over the run lengths and
their starting positions, then round upward to an integer.
These terms all bound mean counts under the uniform message.
Add the three upper counts and clip at binom(2000,v), the
whole shell. Denote this integer by U_v. This upper bound
allows repeated counting when several reasons remove a pair.

At cutoff t<=230 let V_t be the volume of the 2000-bit ball,
and A_t=V_t-sum_(v<=t) U_v. At least Q*A_t ordered pairs of
true and received words remain, for Q=2^1000. For each
received word y write m_y for its number of remaining true
words. Since any two distinct remaining words differ in
every field, the upper bound on the sum of m_y*(m_y-1) is
V_t*D_t from the full-change calculation. These pairs are
all admitted by the actual decoder. Its fallback may increase
the number F_t of received words for which it finds a codeword
within distance t.

For each integer k>=1,

    k*(k+1)*1[m_y>0] >= 2*k*m_y-m_y*(m_y-1).

At a positive integer m_y the difference is
(m_y-k)*(m_y-k-1), which is nonnegative; at zero both sides
vanish. Summing and dividing by Q gives

    F_t/Q >= (2*k*A_t-V_t*D_t/Q)/(k*(k+1)).

For A_t>0, choose k=1+floor(V_t*D_t/(Q*A_t)). In fact the
difference from k to k+1 has the sign of V_t*D_t/Q-k*A_t.
The calculation floors this rational lower bound, clips at
zero and takes its running maximum over the cutoffs, which
is valid because F_t is nondecreasing. Use zero when A_t is
zero. Let the resulting lower integer at cutoff t be G_t
and G_(-1)=0.

The decoder returns a codeword of smallest distance in its
list including fallback. Each received word whose decoded
codeword is at distance t has success contribution
9^(2000-t)/(Q*10^2000). Summation by parts with the
nonincreasing weights 9^(2000-t), and ignoring successful
decodes beyond the cutoff, gives success at least

    sum_(t=0..230) (G_t-G_(t-1))*9^(2000-t)/10^2000.

Part 3 checks the claimed error bound from this exact lower
success bound and the validated witness counts.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The four received fields are copied to padded little-endian
eight-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed data include reduction and square tables, the Golay
and Hamming mapping and candidate tables, and six additional
quadratic-residue tables. The ordinary candidate record has
the systematic error in its low twelve bits, distance in the
next five and cost above those. The additional tables are:

- QR_PARITY_BASE, one parity word per 24-bit systematic input.
- QR_ORBIT_BASE, one word per 24-bit syndrome, giving the
  orbit index in the low nine bits and the permutation index above.
- QR_PERM_BASE, a word at offset k*(6*256)+j*256+v giving
  the image of byte v in byte position j of a noise word
  under permutation k. The RAM representation puts all 24 parity
  bits first and the 24 systematic bits next. In this order, indices
  0..22 retain their raw index, index 23 corresponds to raw 47,
  and indices 24..47 correspond to raw 23..46. To define the
  table apply the raw-coordinate permutation between these two
  reindexing maps. Nonzero image bits occupy disjoint positions
  across input bytes.
- QR_ROW_BASE, one row of 2^23 words per representative,
  each holding two sorted 30-bit entries described above.
- QR_REP_BASE, one canonical representative syndrome per orbit.
- QR_COST_BASE, 49 distance-indexed costs for each of 14 profiles,
  in the same profile order as PROFILE_BASE[24].
- BUDGET_BASE, one word per anchor depth 0 through 100 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (widths 12 and 4) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*(5601)+r for every type p and state r.

Construction.py, qr.py and dynamics.py specify every fixed datum. The
exclusive end of RAM is MEMORY_END in construction.py, within the model's
2^32-word limit. All table addresses are distinct and the
program's temporary words precede these tables.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All 23 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use eight accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=8 words and h=52 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^500. This fits in one word since the
largest possible degree is 24. To multiply a field
value in registers R0..R[k-1] by X^s, emit

    c := shr(R[k-1],h-s)
    for j=k-1 down to 1:
        carry := shr(R[j-1],64-s)
        R[j] := shl(R[j],s)
        R[j] := or(R[j],carry)
    R0 := shl(R0,s)
    R[k-1] := and(R[k-1],2^h-1)
    adr := add(REDUCE_s,c)
    term := load(adr)
    R0 := xor(R0,term)

The removed top coefficients are exactly c; replacing
c*X^500 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=63 byte slots, emit a
load of its source word, a right shift to its byte,
AND with 255, multiply by k, and add the literal
base of the square table for that slot. The table
entry is the k words of (d*X^(8j))^2 for byte d in
slot j. Add the literal word displacement, load
and XOR into the appropriate accumulator for
each of its k words. Store the result into SQ.
By binary linearity this computes the square, at
cost S=b0*(5+3k)+k.

To complete the cube, take the product of SQ and
the same literal source array. Load the k SQ words,
store zero and SQ at multiplication rows zero and
one, and generate rows 2,4,8 by applying SHIFT_1
and storing each result. For every other nonzero
row d, in increasing order, load the row for its
least set bit and the row for the remaining bits,
XOR and store for each word. These are eleven rows,
at four operations per word each. Their contents
are the products of SQ by the binary polynomial d.
Move R registers to zero. For each of n0=125 nibble
slots of the source, in decreasing order, load its
word, shift right, AND with 15, multiply by k and
add the multiplication base. Apply SHIFT_4 followed
by adding the indicated row to R with one address
add, load and XOR per word. Horner's rule yields
the field product. Store R to the designated cube
array. The multiplication cost is

    Mmul=3k+3*(3k+3+k)+44k+n0*(5+3k+3+3k)+k.

Both parts together have cost CubeOps=S+Mmul and give
the exact field cube with zero high padding. The
high padding of the source supplies zeros in the
last partial byte or nibble.

### Bit pieces and initial state

To extract a literal bit piece of width 1 through 64, load its
starting word, and right shift by the within-word offset if
nonzero. For a crossing, also load the next word, shift left
by the complementary offset, and OR. If the part is shorter
than 64 bits, mask it. This costs at most six operations.
Use this for the 32 padded words of Y and P, then store each.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 96. Store zero in the 16
padded incumbent words. Copy Y0 to X0 in sixteen primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear eight
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the eight
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=500, and Depth is a persistent register.

    ANCHOR_VISIT:
        visit the anchor and return to ANCHOR_EXTEND as below
    ANCHOR_EXTEND:
        stop := eq(Depth,R)
        branch if stop to ANCHOR_BACK
        addr := add(NEXT,Depth)
        i := load(addr)
        stop := eq(i,M)
        branch if stop to ANCHOR_BACK
        j := add(i,1)
        store(addr,j)
        Depth := add(Depth,1)
        addr := add(NEXT,Depth)
        store(addr,j)
        addr := add(LAST,Depth)
        store(addr,i)
        FLIP(i)
        jump to ANCHOR_VISIT
    ANCHOR_BACK:
        done := eq(Depth,0)
        branch if done to OUTPUT
        addr := add(LAST,Depth)
        i := load(addr)
        Depth := sub(Depth,1)
        FLIP(i)
        jump to ANCHOR_EXTEND

For FLIP, shift i right by 6 and AND with 63 to obtain
word and bit indices. Shift 1 left by the bit index, add
the literal X0 base to the word index, load, XOR and
store. This costs seven. Each visited subset advances
NEXT before a child, whose own NEXT starts at the next
larger index; LAST records the undone bit. There are
N-1 successful extensions, N failed extensions, N-1
nonroot returns and one root return for
N=sum_(i<=R) binom(500,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base (QR_PARITY_BASE at width 24),
load the parity image, and pack
these images into a eight-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining by adding Depth to the literal BUDGET_BASE and loading
that word. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 560 further primitives:
the 23 cells use at most 23*(6+2+4), and finishing AP
uses 24, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 23 cells in order.
Its frame word slots at offsets 0..5 are Base, U, Rem,
Dist, Rank and Perm. The literal SELECTED slot holds
this cell's selected systematic value. Parent continuation
means the previous block's CELL_NEXT, or ANCHOR_EXTEND
for the first block.

For width m=4 or 12 use the following emitted program.

    CELL_INIT:
        u := extract the received systematic slice
        v := extract the received parity slice
        f := extract the predecessor cube slice
        adr := add(PARITY_BASE[m],u)
        parity := load(adr)
        s := xor(v,f)
        s := xor(s,parity)
        row := shl(s,m)
        base := add(LIST_BASE[m],row)
        store(Base,base)
        store(U,u)
        idx := add(PROFILE_BASE[m],s)
        p := load(idx)
        j := mul(p,5601)
        j := add(j,Remaining)
        idx := add(CAP_OFFSETS[DepthCell-1],j)
        effective := load(idx)
        store(Rem,effective)
        store(Dist,Distance)
        store(Rank,0)
        fall through to CELL_NEXT

    CELL_NEXT:
        rank := load(Rank)
        stop := eq(rank,2^m)
        branch if stop to parent continuation
        addr := add(load(Base),rank)
        record := load(addr)
        cc := shr(record,17)
        needed := add(cc,23-DepthCell)
        current := load(Rem)
        stop := ltu(current,needed)
        branch if stop to parent continuation
        next := add(rank,1)
        store(Rank,next)
        Remaining := sub(current,cc)
        u := load(U)
        emask := and(record,4095)
        selected := xor(u,emask)
        store(SELECTED slot,selected)
        part := shr(record,12)
        d := and(part,31)
        previous := load(Dist)
        Distance := add(previous,d)
        run FINAL if this is the literal final cell
        jump to the next CELL_INIT or, at the final cell,
            back to its CELL_NEXT

Each nested load or arithmetic expression is a separate primitive.
Initialization uses at most 36. The unique stopping attempt
at a visited frame costs at most 11. Each successful attempt
costs at most 28 including its jump, plus FINAL where applicable.

At a width-24 cell change the initialization to use the same three
extracts, the QR_PARITY_BASE lookup and the same two XORs to form
s. Then emit

    addr := add(QR_ORBIT_BASE,s)
    label := load(addr)
    orbit := and(label,511)
    perm := shr(label,9)
    row := mul(orbit,2^23)
    base := add(QR_ROW_BASE,row)
    idx := add(QR_REP_BASE,orbit)
    rep := load(idx)
    packedbase := or(base,shl(rep,32))
    pos := mul(perm,6*256)
    pb := add(QR_PERM_BASE,pos)
    idx := add(PROFILE_BASE[24],orbit)
    p := load(idx)
    cb := add(QR_COST_BASE,mul(p,49))
    packedperm := or(pb,shl(cb,32))
    store(Base,packedbase)
    store(Perm,packedperm)
    store(U,u)
    j := mul(p,5601)
    j := add(j,Remaining)
    idx := add(CAP_OFFSETS[DepthCell-1],j)
    effective := load(idx)
    store(Rem,effective)
    store(Dist,Distance)
    store(Rank,0)
    fall through to CELL_NEXT

At width 24, CELL_NEXT loads rank and returns if rank is at
least 2^24. Otherwise it obtains the next packed entry and
its cost with these instructions.

    framebase := load(Base)
    base := and(framebase,2^32-1)
    addr := add(base,shr(rank,1))
    record := load(addr)
    lane := mul(and(rank,1),30)
    entry := and(shr(record,lane),2^30-1)
    d := shr(entry,24)
    frameperm := load(Perm)
    cb := shr(frameperm,32)
    cc := load(add(cb,d))

Apply the same reserve test using cc, then increment rank
and set Remaining as before. The address halves and the
representative fit in unsigned 32-bit slices. Reconstruct
and permute the raw 48-bit pattern at a successful attempt.

    input := and(entry,2^24-1)
    rep := shr(framebase,32)
    pe := xor(load(add(QR_PARITY_BASE,input)),rep)
    raw := or(pe,shl(input,24))
    pb := and(frameperm,2^32-1)
    move err to zero
    for j=0..5, expanded at emission:
        val := shr(raw,8*j)
        val := and(val,255)
        addr := add(pb,256*j)
        addr := add(addr,val)
        image := load(addr)
        err := xor(err,image)
    e := shr(err,24)
    u := load(U)
    selected := xor(u,e)
    store(SELECTED slot,selected)
    previous := load(Dist)
    Distance := add(previous,d)
    jump to the next CELL_INIT

This applies exactly the recorded raw coordinate permutation
and extracts its systematic part. It charges no more than six
per byte, plus the fixed distance, frame and field updates.
The QR initialization uses at most 60 primitives, including
its three initial extracts, and its single failed attempt at
most 30. Charge at most 90 for these together and at most
100 for each successful attempt including its jump. There is
no QR frame at the final literal cell. All stages use the
six frame words and fewer than 20 scratch registers alongside
the persistent registers specified above.

When an accepted child returns, the parent still has its
Base, U, saved budget and distance and next rank in RAM.
The frame stops when it runs out of entries or finds the
first over-budget record. The costs are sorted, so further
entries would also fail. Each full tuple in the specified
list is therefore visited exactly once.

At FINAL read the 23 selected cell slices and pack them into
the eight padded words of X1, at cost at most
23*(1+4)+8=123. Compute C1=X1^3. Starting from Distance,
for each of the eight words load AP and C1, XOR, popcount
and add to obtain TotalD. The prior accumulated Distance is
the anchor systematic error plus both observations in the
child cells. The residual AP+C1 is P0+L(X0)+C1, so TotalD
is the true full Hamming distance. Compare TotalD with BestD
and branch around replacement unless it is strictly smaller.
If replacing, set BestD := TotalD and copy all sixteen
padded systematic words of X to B. This is at most 320
primitives plus the cube. It retains the earlier candidate
on equality and writes the 1000 systematic bits as prescribed.

## Exact worst-case operation bound

Number the 23 literal cells by i=1..23. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 90 at width 24 and 47 otherwise, and g_i be 100 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
T_i is CubeOps+320 at i=23.
It is zero at all other indices.

For all possible states r=0..5600 define in decreasing cell order

    F_24(r)=0
    U_(i,p)(s)=sum_(c<=s-(23-i))
                 A_(p,c) * (g_i + T_i + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,5600]: U_(i,p)(s)<=F_i(r)-h_i}.

The nonnegative recurrence is nondecreasing in r by induction
and U is nondecreasing in s. The maximum defining E exists
because r itself is allowed. Its search returns the rightmost
valid index, including equality. These are the fixed CAP table entries.
For any actual state r and profile p, the accepted attempts
are exactly the rows up through their recorded cost threshold
E_(i,p)(r)-(23-i). Their continuation states E_(i,p)(r)-c
belong to [0,5600] because costs are positive and successful
attempts satisfy the reserve. At a successful attempt the
remaining program cost is bounded by g_i+T_i+F_(i+1) at that state.
Consequently the actual entire cell traversal costs at most
h_i+U_(i,p)(E_(i,p)(r))<=F_i(r). This proves the bound by
downward induction without distributional hypotheses.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 560+CubeOps. Hence a worst-case
all-input integer bound is

    W = 1000 + sum_(d=0..100) binom(500,d) *
                 (44 + 560 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.

## Scores

The exact counted primitive bound is W=1656702475490230377816015563799334119084235453203098065795607233162982978854442366940820205419897202428804355784407339361400, with
W^500 <= 2^204663. The verifier checks this
inequality and the exclusive RAM address limit.

The exact success calculation above gives error at most
B=11035050186534947/200000000000000000. The declared experiment checks that rational inequality
using integer or fraction arithmetic throughout.
