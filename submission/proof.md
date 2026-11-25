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

The decoder enumerates all anchors x0 of distance d<=113
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
scale S=2^80. Round each local mass down to
floor(S*m_(i,p,c,w)). The empty suffix polynomial has
coefficient S at degree zero. At a step multiply a lower
suffix coefficient by the appropriate local mass integer
and sum contributions, then floor the resulting coefficient
divided by S. Nonnegativity gives, by induction, lower
bounds on S times all the exact P coefficients. Only degrees
through 230 are required, because the scoring contribution
in this proof only uses total noise weights up to 230.

The verifier stores each polynomial with base 2^168, using
one digit per weight. Multiplication by the local mass
integer and a shift by 168*w combines the two factors in
the suffix recurrence. Before dividing at a stage, each
coefficient is at most S^2, since the total local mass is
at most one and every lower suffix coefficient is at most
S. All contributions are nonnegative. There are therefore
no carries among 168-bit digits. Discard digits above the
cutoff, shift the whole integer right by 80, and mask the
low 81 bits of each digit at its original position. This
equals flooring each original digit divided by S: displaced
bits from a following digit begin at bit 168-80=88 and
are erased by the mask. This implements the rounded
probability recurrence using only exact integer operations.

For a true anchor depth d within the radius, start the suffix
at C(d). Let l_(d,z) be its lower integer probability numerator
at child weight z, with denominator S=2^80.
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

The checker evaluates this signed expression exactly for every
mask weight. Put U_w(X)=(1-X)^w*(1+X)^(K-w), and let
P_w(X) be U_w(X) times the reversed intersection polynomial
whose coefficient at degree 2t-d is J_t(d). The coefficient
of P_w at degree 2t-d_s equals the signed expression for
d_s. The first row P_0 is a nonnegative product with
(1+X)^K. The coefficients in the intersection polynomial
are at most V_t and the coefficients of (1+X)^K sum to
2^K. The packed digit stride is chosen beyond V_t*2^K,
so all needed coefficients are read without carry.
Truncating at degree 2t preserves these coefficients.

The polynomial identity (1+X)*P_(w+1)=(1-X)*P_w gives

    P_(w+1)[p] = P_w[p] - P_w[p-1] - P_(w+1)[p-1].

The negative index is zero. Thus the checker propagates
all signed rows in order with exact integer addition and
subtraction. It takes the absolute value only after summing
the signed terms at the desired degree. This proves the
exact weight-specific magnitude used in each rank.

Use the same rank argument as for the per-noise bound.
For one nonzero mask block at weight w, its multiplicity
among first directions is at most 3*binom(M,w). For two
selected blocks at combined weight w the capacity is at
most 9*(binom(K,w)-2*binom(M,w)), since both masks are
nonzero. Keep the capacities at each exact mask weight. When only one
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
on these disjoint sets. Apply the following auxiliary removal rule. A single fair
coin chooses ascending or descending lexicographic message order.
Fix its result across all pairs in the analysis. Remove a true pair
for a run if its single-run competitor is strictly closer to the
received word, or ties the distance and precedes the true message
in the chosen order. This order is used only in the certificate.

For either resulting order, two retained words at the same
received word cannot differ on a proper subset. If there is a
distance decrease over all changed runs, some run decreases it.
If the total change is zero but some run changes the distance,
a run has a negative change. If all runs tie, the one containing
the first differing message bit gives the same lexicographic
preference as the whole change. In this last case the later
word in the chosen order is removed. These arguments apply
also when there is only one changed run.

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
pattern fixes the direction on the run. Let the local true noise weight on the affected coordinates be
w. For fixed actual noise and transmitted message, averaging over
the auxiliary coin weights a candidate of local distance d by
1 if d<w, by 1/2 if d=w, and by 0 otherwise. For any one of the
specified r, r^2 <= 1/2. The three weights are therefore bounded
above by (1/2)*r^(2*(d-w)) for every integer d.

Apply the hyperplane point-mass bound to each actual run direction
and then enlarge to all indexed local error patterns. For each
proposed systematic error the direction is unique. Any additional
zero directions in the enlarged sum only add nonnegative terms.
With the same paired-cell polynomials this bounds the expected
number of preferred competitors over both the uniform message
and the auxiliary coin, conditional on each fixed true noise, by

    ((2/q)^a / 2) * r^(-2*w) * (1+r^2)^u * product H_s(r^2).

The chance of removal is at most the lesser of 1 and that
expected count. Apply min(1,x) <= sqrt(x) for nonnegative x.
The resulting conditional bound is

    2^(-(((500-1)*a+1)/2)) * r^(-w) * (1+r^2)^(u/2)
        * product_(23 pair cells) sqrt(H_s(r^2)).

Here a=1, so the power factor is exactly 2^-250.
Use the rational r=42/100.
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
of removed patterns per true message, averaging also over the
auxiliary orientation coin. Remove a pair of true and received
words if the anchor or child admission test fails or if some
single-run competitor is preferred by the removal rule.
For each realized order, the retained words at any received
word cannot have a proper changed subset between them.

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
These terms bound mean counts under the uniform message
and auxiliary order.
Add the three upper counts and clip at binom(2000,v), the
whole shell. Denote this integer by U_v. This upper bound
allows repeated counting when several reasons remove a pair.

At cutoff t<=230 let V_t be the volume of the 2000-bit ball,
and A_t=V_t-sum_(v<=t) U_v. In expectation over the order,
at least Q*A_t ordered pairs of true and received words remain,
for Q=2^1000. Fix an order. For each
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
vanish. For each order, the number of nonempty retained lists is
at most F_t, which depends only on the actual decoder. Sum the
inequality and then average it over both orders. Dividing by Q gives

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
- BUDGET_BASE, one word per anchor depth 0 through 113 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (widths 12 and 4) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*(5601)+r for every type p and state r.

Construction.py, qr.py and dynamics.py specify the original RAM
regions, including the fixed admission tables. Accelerator.py adds
multiplication-row scratch words starting at their exclusive end,
MEMORY_END. The fixed wide-reduction and wide-square tables follow
that scratch range, ending at FAST_END. Prefix_cube.py places two
state flags, a base cube buffer, two base-tail words, and derivative
lookup work pages at and above FAST_END. Its fixed mode table,
pair of linear maps, unary cube tables and pair contribution
tables follow those temporary words. The formulas below and
the exact address maps in prefix_cube.py specify all entries.
The exclusive end PREFIX_END is 4,287,372,049, at most 2^32.
Fixed and working regions are disjoint.

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
The cube macro uses six vectors of 8 words each and at most six
additional scratch registers. Together with the four persistent
registers this uses at most 58. Other routines' scratch
registers are reused where their lifetimes do not overlap.

### Exact field cube

Write k=8 for the field word count and h=52 for the
number of meaningful high-word bits. Let P be the field modulus
and P_low=P+X^500. The following three stages are the emitted
cube_schedule in accelerator.py, used for a full field cube.
Every cube yields the same field value, so the transmitted code,
candidate admission, costs,
choice of candidate, and output bits in the error proof are
unchanged. The fixed budget and cap tables remain the ones
specified by dynamics.py.

First read the k padded source words into registers. Partition
the source coefficients into successive 21-bit chunks, with a
short final chunk if needed. For each chunk index j and possible
value d, the wide-square table stores k words of (d X^(21*j))^2
modulo P. Extract the chunk with shifts, OR when it crosses a
word boundary, and an AND mask when higher coefficients remain.
Compute its row address using
one multiplication by k and a base add. For each output word
load its table contribution, with one address add if the offset
is nonzero. Move in the first chunk contribution, then XOR the rest. Binary linearity of squaring gives the full
square S in the first monomial vector.

A shift/reduce starts with c=shr(high,h-s) for shift width s.
For each higher word, in decreasing index order, emit a right
shift of its old predecessor, a left shift of the word by s,
and an OR to combine them. Left shift the lowest word by s
and mask the highest word to h bits. The removed coefficients
form c, so the exact reduction is the binary polynomial product
c*P_low. XOR it into the lowest word or words. Its degree is
below the field degree, so no further reduction is needed.

For s=1 the carry c is 0 or 1. An unsigned integer multiplication
by P_low implements its carryless product exactly in one word.

For s=12, REDUCE_WIDE contains one word per possible carry,
holding c*P_low in the binary polynomial ring. Its degree is
at most 32. One address add and one load retrieve that word
before XORing it into the low word.

Use 3 positional four-bit multiplication tables. At group g,
generate four register vectors representing S times X^(4*g+i)
for i=0..3. Each successive vector is obtained from the previous
one by moving it into its named registers and applying the
one-bit shift/reduce. The first vector is the square already
calculated. The zero row of each scratch table is initially
zero and is never written. Traverse the 15 nonzero four-bit
keys in binary reflected Gray order. Start by moving the
first monomial vector into the running row and storing its
k words. Each later key toggles one monomial vector, which
is XORed into the running row before the k words are stored. The index and
toggled coordinate of every Gray step are compile-time constants.
After each group, these stores define all 16 multiples of
S by the associated shifted four-bit polynomial. The physical
nonzero scratch rows are overwritten before each product stage.

Partition the original source into 12-bit windows and
visit them from high to low. Initialize the product registers with the first lookup. Except
at the first window, multiply the product by
X^12 with the wide shift/reduce. For each four-bit part
of the window, extract it and load and XOR the corresponding
positional multiple. The highest partial nibble has only the remaining source
coefficients because the input word is padded. Nonexistent
nibbles are omitted. Horner's rule computes
the product of the source and its square in the field. Store its
k words to the same cube destination as in the original macro.

The expanded instructions are loads, stores and the scalar
operations listed in the model. All extraction indices, table
group bases and word offsets are literal. Dynamic lookup
addresses use the listed multiply/add instructions, remain
inside their tables, and are counted. The emitted sequence has
no data-dependent branches or iteration counts, so counting its
primitive instructions gives CubeOps=6108 for all sources.
Free register moves are excluded under the model. The experiment
checks the exact opcode counts and replays eight representative
field values on the unchanged organizer counter. The replay
shows consistency of the fixed schedule and field arithmetic.
The identity of the cube for every field value follows from
the linear square, modular shift and Horner arguments above.

### Cube reuse within a fixed prefix

The first 17 child cells of field 1 are width-24 cells and
fix its low 408 coefficients. The other six cells determine
the high 92 coefficients. When using the cache for such a
prefix, let a be the first completed field and store its full
cube in BASE_CUBE. Also save the last two words of a in
BASE_TAIL. These include all 92 suffix coefficients. For any
later completed candidate x in this same prefix, d=x XOR a is
supported on coefficients 408 through 499. In the field,

    x^3 = a^3 XOR d^3 XOR a^2*d XOR a*d^2.

The latter two terms form a binary linear map in d for fixed
a. At setup, partition a into literal 17-bit chunks, with a
short last chunk, and the fixed linear table. A row for a
chunk at offset j gives the two field vectors

    (chunk*X^j)^2 * X^408,  (chunk*X^j) * X^816.

XOR the row contributions into two vectors V and U, respectively.
By characteristic two and distributivity, the resulting vectors
are a^2*X^408 and a*X^816. For tail coefficient j from 0 to
91, the derivative basis vector is V XOR U. Advance V by
multiplying by X and U by multiplying by X^2 before the next
coefficient. The same shift/reduce macro is valid at width
two, using the first four carry entries of REDUCE_WIDE.
Thus the vectors at the start of step j are
a^2*X^(408+j) and a*X^(816+2*j). Split these basis vectors
into groups of three, with a shorter last group. The
literal Gray template writes each possible nonzero row as
the XOR of the selected basis vectors. The zero rows start
zero and are never written. The dynamic table pages are
fully regenerated for each new a. Together with saving
the base cube and tail this setup costs 10273 primitives.

To evaluate a cached completion, XOR the final two source
words against the saved words. The known common low prefix
makes this exactly the field difference d, and the extraction
starts at its first possible nonzero coordinate. Load the
saved base cube and XOR the derivative-page contribution for
each three-bit group. This adds the two cross terms in the
displayed identity. To add d^3, partition d into fixed
chunks of widths 10, 10, and eight widths of 9, in order.
Let d_i be each chunk at its literal field coefficient
offset, with zeros elsewhere. Then

    d^3 = XOR_i d_i^3 XOR XOR_(i<j) (d_i^2*d_j XOR d_i*d_j^2).

The unary and pair tables store exactly these fixed functions.
A pair is indexed by the first chunk value in the low bits
and the second in the higher bits. Each XOR loads the eight
consecutive words at its computed row address. All dynamic
indices are bounded by the current chunk widths, so their
computed addresses are within their fixed table intervals.
The straight-line cached completion costs 2348 primitives,
including writing the resulting cube to C1. The static
paths in prefix_cube.py use fewer than 60 registers including
the four persistent ones.

At the literal entry to cell 18, before CELL_INIT, emit

    adr := add(MODE_BASE,Remaining)
    mode := load(adr)
    store(MODE_STATE,mode)
    store(READY_STATE,0)

This executes once per fixed prefix, before any completion.
A return to the same cell's CELL_NEXT does not reinitialize
these state words. MODE_BASE supplies one Boolean for every
incoming state 0 through 5600. The exact work proof defines
the Boolean. It is a fixed table value and controls only the
cube implementation. At each terminal completion, after the
original packing of X1 into its eight words, replace the
ordinary cube call by

    mode := load(MODE_STATE)
    branch if mode to CACHE
    emit the full cube X1 into C1
    jump to SCORE
    CACHE:
        ready := load(READY_STATE)
        branch if ready to REUSE
        emit the full cube X1 into C1
        emit cache_setup from X1 and C1
        store(READY_STATE,1)
        jump to SCORE
    REUSE:
        emit cache_reuse from X1 into C1
        jump to SCORE
    SCORE:
        continue the original terminal distance evaluation

These branch tests use Boolean flags. All cache paths are
within one fixed low prefix until the next cell-18 entry,
so the reuse input conditions hold irrespective of the
order of completed tail words. Both possible modes
produce the original exact cube on each visited word.
They cannot change the candidate list, distances, tie order
or decoded message. The verifier replays cache setup and
twenty reuse completions using the protected counter and
checks the expanded primitive counts and table ranges. The
identities above prove correctness for every field value.

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
23*(1+4)+8=123. Compute C1=X1^3 using the mode program
above. Starting from Distance,
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
To reproduce the fixed admission caps, use the original cube
cost parameter CubeOld=9324 in the cap definition.
Let T_i^0 be CubeOld+320 at i=23 and zero elsewhere.

The unchanged cap recurrence has states r=0..5600.

    F_{24}(r)=0
    U_(i,p)(s)=sum_(c<=s-(23-i))
                  A_(p,c) * (g_i + T_i^0 + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,5600]: U_(i,p)(s)<=F_i(r)-h_i}.

All terms are nonnegative, so F and U are nondecreasing
by downward induction. Each maximum exists since s=r is
eligible. Dynamics.py derives these exact cap values, using
the original cube cost solely for this data definition.
This is also the cap table used by the true-path probability
calculation. No cap is changed by accelerating the cube.

For the actual operation bound, first compute two arrays for
each of the suffix cells i=18..23. Their common boundary is
G^plain_24(r)=G^reuse_24(r)=0. In each mode define

    G_i^mode(r) = h_i + max_p sum_(c<=E_(i,p)(r)-(23-i))
                    A_(p,c) * (g_i + T_i^mode
                                  + G_(i+1)^mode(E_(i,p)(r)-c)).

For both modes T_i is zero for i<23. For i=23 use
T^plain_23=320+6108+3 and T^reuse_23=320+2348+6.
The extra three in the plain case cover its flag load,
branch and jump. In cache mode, every completion costs
at most six flag and control primitives, and a cached
completion uses 2348 cube primitives. For the first
completion charge an additional P=6108+10273 once per
fixed prefix for its full cube and cache setup. This
bound charges the cached completion work even on the
first completion, so it is an upper bound there as
well. Charging P once is valid even for prefixes with
no completed candidate.

Accordingly MODE_BASE at state r stores one exactly
when G^reuse_18(r)+P < G^plain_18(r). The literal
entry's four setup instructions give the bound

    G_18(r)=4+min(G^plain_18(r), G^reuse_18(r)+P).

For i=17 down to 1, use the same recurrence with
G_(i+1) and T_i=0. The exact admission cap remains
E_(i,p)(r), and each continuation state after a
successful attempt stays in 0..5600 by the reserve
and positive cost. For each profile the sum bounds
the attempted row costs and recursive work. The
existing h_i covers initialization and the final
stopping attempt at every frame. This proves the
all-input bound by backward induction with a distinct
cost for each cache mode. No branch count relies on
how often that mode occurs under the channel.

The anchor traversal, parity setup and full cube
charge remain at most 44+560+6108 for each visited
anchor, with initialization and output bounded by
1000. Thus the exact primitive bound is

    W = 1000 + sum_(d=0..113) binom(500,d) *
                 (44+560+6108+G_1(C(d))).

The declared exact recurrence regenerates the inherited
caps and the new mode table, calculates W, and checks
the exact exponent inequality.

## Scores

The exact counted primitive bound is W=12090556300019517937004969428759815635665987178972565222314628033298170323706619569510172111815128586782289195893137422017941174208690, with
W^9677 <= 2^4278108. The verifier checks this
inequality and the exclusive RAM address limit.

The exact success calculation above gives error at most
B=13793812699701584233989266086719410691211/250000000000000000000000000000000000000000. The declared experiment checks that rational inequality
using integer or fraction arithmetic throughout.
