# Ten degree-100 cyclic fields with exact weight-indexed admission

## Encoder and decoder

Use the polynomial basis of F2[X]/(P), where
P=X^100+X^6+X^5+X^2+1. The declared field check applies
the Rabin criterion with prime divisors 2 and 5 of 100.
Split the message into ten 100-bit words b0 through b9.
Transmit these ten words, followed by the parity fields
p_i=L(b_i)+b_(i-1)^3 with cyclic indices. The 1000-bit
systematic prefix gives injectivity and the twenty fields give
2000 channel bits.

In each field L is block diagonal with four width-24 cells and
one width-4 cell in increasing bit order. For a width-24
input u, the low 23 parity bits are the binary polynomial
remainder of u*X^23 modulo the degree-23 G=9205487 (in
low-coefficient-first order). The high parity bit is the parity
of the sum of Hamming weights of u and the low parity bits.
For each 4-bit Hamming cell the parity columns are 14,13,11,7.
The five cell maps are linear. Construction.py fixes their integer layout.

The decoder enumerates all anchors x0 of distance d<=100
from received Y0. For each use entry C(d) at index d of construction.py BUDGETS as the initial
remaining-cost state. Visit the child fields j=1 through 9 in order. In each field,
form the local systematic and adjusted parity observations u,v
from received Yj and Pj+x_(j-1)^3. At width m define the local
syndrome s=v+L_m(u). For each error e on the systematic
part, the paired distance and unnormalized weight are

    d_s(e)=wt(e)+wt(L_m(e)+s),
    A_s(e)=9^(2m-d_s(e)),   S_s=sum_e A_s(e).

Assign cost c_s(e)=ceil(32 log2(S_s/A_s(e))). The preloaded
integer cost is bit_length((S_s^32-1)//A_s(e)^32), since
the quotient is below 2^c exactly when S_s^32<=2^c*A_s(e)^32.
Costs are positive because every complete row has more than
one positive term.

All candidates of width 4 appear in fixed rows
ordered by (distance, numerical systematic error). In width 24
use a row of the first 4194304 noise patterns, through an explicit
syndrome orbit map and coordinate permutation defined below.
This preserves the distances and the syndrome. The selected
cell is u+e, where e is the systematic part of its noise pattern
in original coordinates. All rows are in nondecreasing distance
and thus nondecreasing cost.

Depth-first search the 45 child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(45-i), where r is
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

For each representative, the fixed row lists all raw 48-bit
noise patterns with that syndrome, in increasing (distance,
numeric raw noise), keeping only the first 4194304. For each of
the 2^24 systematic errors its parity error in the row is
L_m(e)+r, so the full list and ordering are explicit. Its
normalizer S_r sums over the entire coset before truncation.
Each retained record is the raw 48-bit error in low bits and
its cost in upper bits. The cost fits in 16 bits, since
A>=1 and S_r<=2^24*9^48. The runtime uses the stored
permutation to map a representative's error back to syndrome s.
A coordinate permutation preserves weight and is a bijection,
so for any weight with coset frequency h_s(d) the number of
retained records of that weight is the smaller of its frequency
and the unfilled prefix capacity at lower weights. This is true
also when the row ends inside a tie at a fixed weight.

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
and the anchor. Each width-m cell has b_i=2m observed bits. Let H_(i,p,c,w) be the number of retained true noise
patterns of weight w and nominal cost c across all
syndromes of its profile type p. Multiply the retained
count for one syndrome by its profile multiplicity.
The QR orbit permutation preserves total weight and cost,
even if truncation cuts the last retained weight tie.
The channel mass for a term (p,c,w) is

    m_(i,p,c,w) = H_(i,p,c,w)*9^(b_i-w)/10^b_i.

For each i, use the same effective cap E_(i,p)(r) and
remaining state r as the decoder. With s=E_(i,p)(r),
a true pattern of cost c passes precisely when it is retained in
its row and c<=s-(N-i), where N is the fixed number of
child cells in the program. The next state is s-c.
On the true path the predecessor adjustment leaves the actual
local noise syndrome. Define
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

## Competitors changing all ten fields

Write q=2^100 and Q=q^10. For nonzero e, the cube
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
annihilate trace on a spanning set. The equation for a nonzero multiplier e^-3 has degree three after inversion, so each nonzero normal mask determines at most three directions e.

A full-change competitor is specified by nonzero
directions e_i in all ten fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the ten translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 1000
parity bits can be evaluated by expanding the
ten affine indicators. Multiplication by their
common point mass 1024/Q leaves 1/Q times a sum
of 1024 terms. The empty mask term is the
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

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(100*h,w).

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
most binom(100*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(100*(10-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(100*(10-k),b) *
    H_k(t-a-b, binom(100*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(10,k) selections of masks, and add
the 2000-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_fixed(t). The expected
number of full-change competitors within total
t for each fixed true noise is at most D_fixed(t)/Q.
The checker verifies D_fixed(t) is nondecreasing through t=230.

## Averaged full-change ball intersections

The occupancy calculation uses the total number of ordered
full-change pairs with a common received word within distance
t of both codewords. For cutoffs 217 through 226 we bound this
count directly by summing ball intersections over all actual
true noise patterns in the ball.

Let M=100, G=10, K=1000 and Q=2^K. Write J_t(d) for the
intersection size of two radius-t balls at distance d.
On the d differing coordinates choose a differences from
the first word and d-a from the second. The remaining equal
coordinates can contribute at most t-max(a,d-a) differences.
Therefore

    J_t(d)=sum_(a=0..d) binom(d,a) *
                 sum_(b=0..t-max(a,d-a)) binom(2000-d,b).

A negative inner limit means zero. Thus J_t(d)=0 for d>2t
and 0<=J_t(d)<=V_t, the volume of one ball.

Fix all G nonzero systematic changes and let d_s be their
combined weight. Their parity difference law averaged over
uniform true fields has uniform density term 2^-K, plus one
signed binary character of that magnitude for each nonempty
choice of the supported parity masks. This follows from the
product of affine difference laws proved above. A selected
mask of total weight w has term, in absolute value before
dividing by 2^K, bounded by

    |sum_(p=0..2t-d_s) J_t(d_s+p)*K_p(w)|,

where K_p(w) is the X^p coefficient of
(1-X)^w*(1+X)^(K-w), obtained by summing the character over
each parity shell. Each fixed translation contributes only a
sign of absolute value one.

For j=min(w,K-w) pair j factors (1-X)(1+X). Taking
absolute values of coefficients bounds the resulting signed
polynomial by

    U_j(X)=(1+X^2)^j*(1+X)^(K-2j).

This is coefficientwise nonincreasing in j, since the
increase substitutes (1+X^2) for (1+X)^2. Round j down to
its bin endpoint in MASK_BINS. The nonnegative convolution
of the U_j coefficients with J_t(d_s+p) bounds the signed
term for every weight in that bin.

Write N_k(d) for the number of d-weight patterns on k
M-bit blocks with each block nonzero. Inclusion-exclusion gives

    N_k(d)=sum_(i=0..k) (-1)^(k-i)*binom(k,i)*binom(i*M,d).

In particular N_0(0)=1 and N_0(d)=0 for d>0. Fix a subset
of k selected parity blocks. Its masks have distinct
predecessor systematic fields. The total weight-w mask
capacity among actual direction tuples is at most
mu^k*N_k(w), where mu=3 in this even degree, by the cubic-normal
bound. Sum those capacities over each bin. For complement
systematic weight b there are N_(G-k)(b) possible nonzero
patterns. The selected fields supply at most N_k(d_s-b)
choices. Sum that many largest upper terms within the
binned mask capacities, multiply by N_(G-k)(b), and sum over
b for which those counts can be positive. Multiply by
binom(G,k) for all chosen subsets. Sum over k=1..G and
d_s=G..2t to obtain the integer correction R_t.

The checker uses a nonnegative packed convolution with
J_t in reverse order, so the coefficient at degree 2t-d_s
is the bound for that total systematic weight. The
nonnegative coefficients of U_j sum to 2^(K-j)<=2^K.
Every digit after multiplying by the reversed J_t
polynomial is at most V_t*2^K. The declared digit stride
has capacity beyond this, so no carries affect any
included coefficient. Truncating at degree 2t preserves
the coefficients used. Bin capacities and rank sums are
computed with exact integers.

The empty mask term is at most V_t^2 before the 2^-K
factor, by including all systematic differences. Each
ordered pair of errors in one ball contributes once at
their XOR difference. Averaging over Q=2^K true words
cancels the parity factor, so the ordered full-change pair
count is at most V_t^2+R_t. It is also at most
V_t*D_fixed(t) by the per-noise proof. At the specified
cutoffs take the lesser of D_fixed(t) and
ceil((V_t^2+R_t)/V_t), and at others use D_fixed(t). Take
the running maximum through t and call it D_t. It still
satisfies the bound V_t*D_t for the total pair count.
Part 1 checks every declared integer in pairs.json.

## Proper-run competitors

A proper nonempty changed set is the disjoint union of its
cyclic runs. Each run of length a changes a systematic fields
and the a+1 parity fields starting with the first changed field.
The affected coordinates of different runs are disjoint.
Order the messages lexicographically in their systematic bit
order. Remove a true pair for a run if its single-run competitor
is strictly closer, or is at equal distance and lexicographically
earlier. If a proper-set competitor is strictly closer, at
least one run is strictly closer. If the proper competitor
ties and is lexicographically earlier, either some run is
strictly closer or each run ties and the one containing the
first changed systematic bit is lexicographically earlier.
Thus two retained words at a given received word cannot have
a proper changed set: remove the farther one or, on a tie,
remove the lexicographically later one. The decoder's
tie order can be any fixed order, because its output
likelihood under the uniform prior depends only on distance.

Fix a run and its true noise. In the first systematic and parity
pair the differences are e and L(e), since its predecessor
is unchanged. The actual syndrome of each cell is invariant.
For fixed syndromes use the local polynomials H_s(X) with
coefficient h_s(d) for the number of possible weight-d
resulting first-pair errors. There are five cell factors, four of width 24 and one of width 4. Each possible
systematic error determines its direction uniquely. The remaining
a-1 systematic and a parity fields have u=(2a-1)*100 bits.
Write M for the field degree, p=a*M, n=(a-1)*M, and u=p+n.
For a fixed nonzero first direction and fixed true noise, use
D_(a,t)/2^(p+1) to bound the expected number of other
directions whose candidate outside-pair distance is below t,
plus half the expected number whose outside distance is t.
This is the average of the ball indicators at t and t-1.

Each outside parity field affected by a nonzero predecessor
has an independent affine hyperplane difference as above.
Expand the product of the p-bit parity restrictions. The
doubled empty-mask contribution is the sum of Hamming ball
volumes at radii t and t-1 on the p+n bits. A signed parity
term at remaining radius l has magnitude
abs(K_l(w)+K_(l-1)(w)), where w is its mask weight and
K_(−1)(w)=0. This is the interval sum for two balls with
the same translation. All terms have the common expected
scale 2^(-p-1).

In a mask term let k be the number of selected masks
belonging to the other a-1 changed directions. There are
binom(a-1,k) choices of that subset. For k>0 their
combined weight-w multiplicity is bounded by

    mu^k * sum_(j=1..k) (-1)^(k-j)*binom(k,j)*binom(j*M,w),

where mu is the bound on preimages of a nonzero normal in
this field, here mu=3. The coefficient comes from requiring all k mask
blocks to be nonzero. If the fixed first mask is absent,
rank these multiplicities by abs(K_l(w)+K_(l-1)(w)).
If the first mask is present, rank by that magnitude
maximized over its possible added weight 1..M.
At candidate systematic weight d on the k masked fields
there are at most binom(k*M,d) choices of their directions.
Sum that many highest magnitudes. On the other systematics
there are at most binom(n-k*M,d') choices of weight d'.
For cutoff t use l=t-d-d' and sum over nonnegative
d,d'. If k=0, the only nonempty mask is the fixed first
one, with magnitude maximized over weights 1..M.
These computations yield D_(a,t) uniformly for every
nonzero fixed first direction and every actual noise.

For each tuple of first-pair cell syndromes let H be the
product of their full weight enumerators. For fixed local
true noise weight w, the number of local true noise patterns
with this tuple is the coefficient of X^w in H(X)*(1+X)^u.
The actual first-pair weights and the u unrestricted bits
are counted jointly here.

On a distance tie, swap the true word and its run competitor
at the same received word. This preserves the run, the local
true total weight and the noise outside the run. It also preserves
each first-pair syndrome, because the first predecessor is
unchanged and the local differences have form e,L(e).
The swap reverses lexicographic order and is a fixed-point-free
involution. In each of these conditioned classes exactly half
the tied ordered pairs are oriented toward the earlier word.
Strict improvements contribute their full counts. The number
of removed pairs is at most the number of such oriented
ordered pairs, and at most the number of true pairs.

For each actual first-pair noise and each nonzero first
direction, the uniform D_(a,t)/2^(p+1) bound above applies
to the outside changes. If H_d is the coefficient of X^d,
at most H_d first directions give a first-pair candidate
weight d. The one unchanged first direction can be omitted;
leaving H_d as the multiplicity only enlarges the upper
bound on the sum over nonzero directions. Therefore, summed
over any group with fixed syndrome tuple and fixed w, the
fraction of true pairs removed is at most

    min(1, 2^(-p-1) * sum_(d<=w) H_d * D_(a,w-d)).

Here the strict and half-tie expectations are uniform over
the actual noise before grouping, and the orientation
identity holds after grouping. Because the conditional
candidate coefficient bounds are nonnegative, for any
z>0 the sum of coefficients is further bounded by

    H(z) * z^(-w) * max_(0<=t<=w) D_(a,t)*z^t.

The checker uses the literal integer radix grids a0
for each run length from RUN_RADICES. At radix a0,
put z=(a0/100)^8. It partitions the first-pair syndrome tuples
by a verified upper bin for H(z). For each component local
polynomial h, the checker evaluates h(z) exactly and rounds it
up to an integer multiple of 2^-76. Let N0 be the product
of the five rounded numerators for the tuple. For

    k = ceil(log2(N0^8)) - 8*76*5

the inequality H(z)<=2^(k/8) follows by comparing eighth
powers; ceil(log2(N0^8)) is (N0^8-1).bit_length().
Every local enumerator is nonzero and z>0, so N0 is positive.

The partial histograms group all syndrome tuples in this bin
without assuming that the same polynomial occurs within one
bin. For repeated widths the checker lists multisets of the
component profile types; it multiplies by the exact number
of orderings and by their syndrome multiplicities. For each
bin it sums H(X) over its tuples and stores the polynomial
A_k(X). This sum is a nonnegative coefficient count of a
subset of all 200-bit true noise patterns. The checker
verifies the summed bin coefficients are binom(200,d).

For each bin k the number of true local noises of weight w
is the coefficient of X^w in A_k(X)*(1+X)^u. Representing
polynomials in base 2^(8*stride) packs each coefficient into
one digit. When forming A_k, the stride has room for an
integer of size 2^200; no intermediate coefficient can exceed
the full 200-bit count. When multiplying by the free
binomial factor each output digit through the cutoff is at
most binom(200+u,w), and the stride is chosen from the largest
of these binomial counts. Reduction modulo the corresponding
power of the base truncates the polynomial at the cutoff.
All smaller coefficients therefore remain exact.

To compute the bound factor, use the integer recurrence

    Z_w = max(D_(a,w)*a0^(8*w), 100^8*Z_(w-1)), Z_(-1)=0.

Then Z_w/a0^(8*w) is exactly
z^(-w)*max_(t<=w) D_(a,t)*z^t. For k=8e+j, 0<=j<8,
round 2^(j/8) upward to denominator 2^64 by an integer
eighth root and exact perfect-power comparison. Multiply by
2^e and use this upper bound on 2^(k/8) in the displayed
proper-run probability. Round the probability upward at
scale 10^45 and clip at one. The clipped upper bound
depends only on (k,w). Multiplying it by the bin's true
noise count is valid after summing over its tuples, which
preserves the tie involution within each. Sum the weighted
counts over bins and take the least of the computed sums
at each w. Each whole sum is an upper bound at that w,
so taking the least is valid.

Convolve this count bound in w with binom(2000-(2a+1)*100,t-w)
for the unaffected noise positions. The factor 9^(2000-t)/10^2000
is the probability of one pattern at true total weight t.
This gives the proper-run removal bound J_(a,t) used below.
All comparisons and arithmetic are exact integers or rationals
with the indicated upward rounding.

The run-removal totals are verified in two portions. For
each of the 231 total noise weights, proper.json stores
the sum of the bounds for lengths 1 through 5 and all
ten start positions, scaled by 10^45. Part 1 verifies
these integers and the signed pair cutoffs. Part 3
recomputes lengths 6 through 9, adds the validated early
numerators and divides by the common scale with upward
rounding. The resulting integer is the proper-run term
in the occupancy bound. All convolutions before this
rounding are exact, and the two ranges cover each proper
run length once.

## Occupancy after removals by exact noise weight

At each total noise weight v, upper-bound the mean number
of removed patterns per true message. Remove a pair of true
and received words if the anchor or child admission test fails,
or if some single-run competitor is closer, or ties and is lexicographically earlier. Two retained
words that differ on a proper subset cannot coexist, since the
farther is removed, or the later one on a tie.

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
The twenty received fields are copied to padded little-endian
two-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed data include reduction and square tables, the Hamming mapping and candidate tables, and four additional
quadratic-residue tables. The ordinary candidate record has
the systematic error in its low twelve bits, distance in the
next five and cost above those. The additional tables are:

- QR_PARITY_BASE, one parity word per 24-bit systematic input.
- QR_ORBIT_BASE, one word per 24-bit syndrome, giving the
  orbit index in the low nine bits and the permutation index above.
- QR_PERM_BASE, a word at offset k*(6*256)+j*256+v giving
  the image of byte v in byte position j of a raw noise word
  under permutation k. Nonzero image bits occupy disjoint
  positions across input bytes.
- QR_ROW_BASE, one sorted 4194304-entry list per representative,
  with each low 48-bit raw noise word and its cost above it.
- BUDGET_BASE, one word per anchor depth 0 through 100 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (width 4) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*(5501)+r for every type p and state r.

Construction.py, qr.py and dynamics.py specify the RAM through
its original exclusive endpoint MEMORY_END with the new row
capacity. Pair_cube.py specifies the 28 fixed pair tables that
start there. Their exclusive ending address PAIR_END fits the
2^32-word limit. Both original and added tables follow every
mutable decoder region.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All 45 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The pair cube uses two input words, seven shifted low chunks,
a shifted high chunk per needed width and index, two output
accumulators and temporary scalars. The declared schedule count
verifies at most 64 names even together with the four persistent
registers. Scratch registers in the rest of the decoder can be
reused at block boundaries.

### Exact degree-100 field cube by pair lookup

Split the operand into eight consecutive chunks of widths
13,13,13,13,12,12,12,12 from the low polynomial coefficients.
Write x_j for chunk j shifted to its original coefficient
positions. Squaring is binary-linear. Distributing the product
of x with its square gives

    x^3 = XOR_j x_j^3 XOR XOR_(i<j)(x_i^2*x_j XOR x_i*x_j^2).

For each ordered index pair i<j preload a contiguous two-word
table. The key is the unsigned i chunk at the bottom and the j
chunk above. Each value contains the corresponding cross term
reduced in the fixed degree-100 field. Include the individual
cube x_i^3 in that table exactly when j is its designated unary
neighbor, and include x_j^3 exactly when i is its neighbor.
Choose neighbor 1 for chunk 0 and neighbor 0 for every other
chunk. This assigns exactly one table to each unary term.
Equivalently, compute an entry as (x_i+x_j)^3 XOR the individual
cubes assigned elsewhere. The exact sum of the 28 table values
therefore equals the cube for every possible source.

The emitted macro loads the two padded source words. Extract
each chunk at its fixed bit offset. For each chunk used as a
low part keep its value shifted left by one. For each high
chunk j and width v of a lower chunk, keep its value shifted
left by v+1. These shifts encode the two-word stride. For pair
i<j, OR the proper two shifted parts and add the pair's literal
base address. Load two words, adding one to the address for the
high half. Move the first pair into the accumulator and XOR
each subsequent pair. Write the two result words. The logical
coefficients occupy at most 100 bits, so the high output word
is correctly padded.

This fixed schedule has 58 loads, 56 adds, 54 XORs, 28 ORs,
17 left shifts, 6 masks, 6 right shifts and 2 stores.
CubeOps=227 is the same on every input. Table access addresses
are within each declared range and the whole table end is
checked. The experiment replays 19 source operands on the
byte-identical reference counter, with exact preload entries.
The maximum cost and field identity hold on every operand
by the literal schedule and the pair expansion.

### Bit pieces and initial state

To extract a literal bit piece of width 1 through 64, load its
starting word, and right shift by the within-word offset if
nonzero. For a crossing, also load the next word, shift left
by the complementary offset, and OR. If the part is shorter
than 64 bits, mask it. This costs at most six operations.
Use this for the 40 padded words of Y and P, then store each.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 120. Store zero in the 20
padded incumbent words. Copy Y0 to X0 in four primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear two
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the two
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=100, and Depth is a persistent register.

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
N=sum_(i<=R) binom(100,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base (QR_PARITY_BASE at width 24),
load the parity image, and pack
these images into a two-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining by adding Depth to the literal BUDGET_BASE and loading
that word. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 69 further primitives:
the five cells use at most 5*(6+2+4), and finishing AP
uses 6, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 45 cells in field and within-field order.
Its frame word slots at offsets 0..5 are Base, U, Rem,
Dist, Rank and Perm. The literal SELECTED slot holds
this cell's selected systematic value. Parent continuation
means the previous block's CELL_NEXT, or ANCHOR_EXTEND
for the first block.

For width m=4 use the following emitted program.

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
        j := mul(p,5501)
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
        needed := add(cc,45-DepthCell)
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
        run FIELD after each of the eight nonfinal fields
        run FINAL if this is the literal final cell
        jump to the next CELL_INIT or, at the final cell,
            back to its CELL_NEXT

Each nested load or arithmetic expression is a separate primitive.
Initialization takes at most 18 for the three extracts and
17 for their following lines. The unique stopping attempt
uses at most 11, giving a total of 46. Each successful attempt
takes 14 through the remaining-state update and 9 through
the jump, for 23 plus FIELD or FINAL where applicable.

At a width-24 cell change the initialization to use the same three
extracts, the QR_PARITY_BASE lookup and the same two XORs to form
s. Then emit

    addr := add(QR_ORBIT_BASE,s)
    label := load(addr)
    orbit := and(label,511)
    perm := shr(label,9)
    row := mul(orbit,4194304)
    base := add(QR_ROW_BASE,row)
    pos := mul(perm,6*256)
    pb := add(QR_PERM_BASE,pos)
    store(Base,base)
    store(Perm,pb)
    store(U,u)
    idx := add(PROFILE_BASE[24],orbit)
    p := load(idx)
    j := mul(p,5501)
    j := add(j,Remaining)
    idx := add(CAP_OFFSETS[DepthCell-1],j)
    effective := load(idx)
    store(Rem,effective)
    store(Dist,Distance)
    store(Rank,0)
    fall through to CELL_NEXT

The two loop tests and ordering in CELL_NEXT are the same, using
rank limit 4194304 and cc=shr(record,48). After the successful
cost test, increment rank and compute Remaining as before.
The remaining action replaces that block's ordinary extraction
and scoring updates with

    raw := and(record,2^48-1)
    d := popcount(raw)
    pb := load(Perm)
    move err to zero
    for j=0..5, expanded at emission:
        val := shr(raw,8*j)
        val := and(val,255)
        addr := add(pb,256*j)
        addr := add(addr,val)
        image := load(addr)
        err := xor(err,image)
    e := shr(err,23)
    e := and(e,2^24-1)
    u := load(U)
    selected := xor(u,e)
    store(SELECTED slot,selected)
    previous := load(Dist)
    Distance := add(previous,d)
    jump to the next CELL_INIT

This applies exactly the recorded raw coordinate permutation
and extracts its systematic part. It charges no more than six
per byte, plus the fixed distance, frame and field updates.
The QR initialization uses at most 18 for three extracts,
4 for parity and syndrome formation, and 20 for the remaining
initialization lines. Its stopping attempt uses at most 11.
Charge 53 for these together. A successful attempt uses 14
through the remaining-state update, 3 before the byte loop,
36 in the byte loop and 8 through the jump, totaling 61. Each field boundary falls at an ordinary width-4 cell. All stages use the
six frame words and fewer than 20 scratch registers alongside
the persistent registers specified above.

When an accepted child returns, the parent still has its
Base, U, saved budget and distance and next rank in RAM.
The frame stops when it runs out of entries or finds the
first over-budget record. The costs are sorted, so further
entries would also fail. Each full tuple in the specified
list is therefore visited exactly once.

At FIELD and FINAL, pack the five selected cell slices from
the current completed field into its two padded Xj words,
with at most 5*(1+4)+2=27 primitives. Compute Cj=Xj^3.
At FIELD this cube supplies the next field; its total routine
charge is at most CubeOps+30 including packing and transfer. At FINAL add the distance of
the anchor parity to Distance by loading AP and C9 round by
round over two words, XORing, popcounting and adding.
The residual AP+C9 is P0+L(X0)+C9. Distance already includes
the anchor systematic distance and all nine visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all twenty systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 84 primitives. Packing costs 27,
the two parity words cost ten, the comparison and branch
cost two, copying the twenty words costs 40, and the
remaining charge covers control transfer.
## Exact worst-case operation bound

Number the 45 literal cells by i=1..45. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i^0 be 53 at width 24 and 46 otherwise, and g_i^0 be 61 at
width 24 and 23 otherwise. These are the actual decoder setup and
row step charges, including one stopping attempt in setup. Use T_i^0
equal to 227+30 at i=5,10,15,20,25,30,35,40, equal to 227+84 at
i=45, and zero elsewhere. Thus the cap generation uses the same
pair-cube and packing charges as the emitted traversal.

The cap recurrence with the new groups has states r=0..5500.

    F_{46}(r)=0
    U_(i,p)(s)=sum_(c<=s-(45-i))
                  A_(p,c) * (g_i^0 + T_i^0 + F_(i+1)(s-c))
    F_i(r)=h_i^0+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,5500]: U_(i,p)(s)<=F_i(r)-h_i^0}.

All terms are nonnegative, so F and U are nondecreasing
by downward induction. Each maximum exists since s=r is
eligible. Dynamics.py derives these exact cap values. The same effective cap
table is used by the true-path probability calculation and in the
actual decoder. The rounded or approximate values used in research
to propose the literal anchor budgets are absent from this recurrence.

For the actual operation bound define G_{46}(r)=0. Let h_i
and g_i be respectively 53 and 61 at width 24, and 46 and 23
otherwise. Let T_i^* be CubeOps+30 at the eight intermediate
field boundaries, CubeOps+84 at i=45, and zero elsewhere.
For each cell i, set

    U^*_(i,p)(s)=sum_(c<=s-(45-i))
                   A_(p,c) * (g_i + T_i^* + G_(i+1)(s-c))
    G_i(r)=h_i+max_p U^*_(i,p)(E_(i,p)(r)).

For any actual state r and profile p, the emitted list traversal
accepts exactly the entries through its fixed cap minus the
reserve. The resulting states E_(i,p)(r)-c stay in [0,5500].
By downward induction the accepted attempts cost at most the
terms in U^*, and h_i covers the setup and the stopping attempt.
Thus G bounds every execution. This direct calculation uses
the fixed admission caps and the actual primitive work, without a
worst-case claim based on inputs sampled from the channel.

The anchor subset traversal uses at most 44 primitives per visit,
and its parity residual setup and cube at most 69+CubeOps.
Together with fixed initialization and output the resulting
all-input bound is

    W = 1000 + sum_(d=0..100) binom(100,d) *
                 (44 + 69 + CubeOps + G_1(C(d))).

The declared experiment regenerates the fixed caps, computes this
recurrence, and checks the exact exponent inequality.

## Scores

The exact counted primitive bound is W=21306706949390502987836000332430252282142597733, with
W^10 <= 2^1539. The verifier checks this
inequality and the exclusive RAM address limit.

The exact success calculation above gives error at most
B=1404719304643696798689007504653/10000000000000000000000000000000. The declared experiment checks that rational inequality
using integer or fraction arithmetic throughout.
