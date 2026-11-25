# Nine degree-111 cyclic fields with punctured anchor parity

## Encoder and decoder

Use the binary polynomial basis of the field defined by
P(X)=X^111+X^11+X^2+X+1. The verifier checks the Frobenius
criterion, including the two prime factors 3 and 37 of 111.
Split the first 999 message bits into nine consecutive 111-bit
field elements x0,...,x8. The last message bit is b. In increasing
bit order divide each x field into four width-24 cells, one
width-12 cell and three width-1 cells. Let L be the blockwise
linear map using the extended quadratic-residue map below at
width 24, the extended Golay map at width 12, and identity at
width 1. For the Golay map take the remainder of u*X^11 modulo
1+X^2+X^4+X^5+X^6+X^10+X^11. The high parity bit makes the
total parity of u and the remainder even.

Before puncturing, the parity field at j is L(xj)+x_(j-1)^3,
with indices modulo nine. Delete the three highest bits of parity
field zero, which are its three width-1 cells. No systematic
coordinate is deleted. The encoded layout is the 1000 original
message bits, then 108 bits of parity field zero and the eight
remaining 111-bit parity fields, then four more copies of b. This
carries all 1000 bits on exactly 2000 coordinates, with independent
BSC observations on the 1995 cyclic coordinates and on b's five
coordinates.

Decode b by majority. For the cyclic part start with a zero
incumbent and enumerate all x0 at systematic distance d<=25
from its received systematic word. Use the integer C(d) stored
in BUDGETS as the initial remaining cost. Traverse the other
eight fields in order, visiting each field's eight cells. A
child's observations are its received systematic slice u and
its parity slice v after subtracting the preceding candidate's
cube. Its syndrome is s=v+L_m(u), for the appropriate local map.
The true noise candidate e pairs systematic error e with parity
error L_m(e)+s. Its local distance is the sum of those two weights,
its candidate systematic slice is u+e, and its integer cost is

    c_s(e) = ceil(16*log2(S_s/A_s(e))),
    A_s(e) = 9^(2m - local distance),
    S_s = sum_e A_s(e).

The exact precomputed cost equals
bit_length((S_s^16-1)//A_s(e)^16). Widths 12 and 1 use complete
rows ordered by distance then systematic error. Width 24 uses
the retained orbit rows described below. The decoder visits the
64 literal cells in a depth-first search. At cell i and current
state r it uses the profile-dependent effective cap E_(i,p)(r).
A listed entry with cost c passes when c<=E_(i,p)(r)-(64-i).
Its next state is E_(i,p)(r)-c. The positive costs justify the
reserve and let a frame stop at its first over-budget record.
For a complete candidate add the missing anchor parity distance
on the 108 surviving parity bits. Replace the incumbent on
strict improvement. Every visited word is a codeword of the
punctured cyclic code. The fallback and fixed tie order define
the decoder for every received word.

## Extended 48-bit cell data and its orbits

The QR generator G is the binary polynomial with coefficient
bit mask 9205487, with the low bit its constant term.
For a 24-bit slice u, take the remainder of u*X^23 modulo G.
The high parity bit makes the total parity of u and the
23-bit remainder even.

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
numeric raw noise), keeping only the first 2^23. For each of
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
scale S=2^48. Round each local mass down to
floor(S*m_(i,p,c,w)). The empty suffix polynomial has
coefficient S at degree zero. At a step multiply a lower
suffix coefficient by the appropriate local mass integer
and sum contributions, then floor the resulting coefficient
divided by S. Nonnegativity gives, by induction, lower
bounds on S times all the exact P coefficients. Only degrees
through 230 are required, because the scoring contribution
in this proof only uses total noise weights up to 230.

The verifier stores each polynomial with base 2^104, using
one digit per weight. Multiplication by the local mass
integer and a shift by 104*w combines the two factors in
the suffix recurrence. Before dividing at a stage, each
coefficient is at most S^2, since the total local mass is
at most one and every lower suffix coefficient is at most
S. All contributions are nonnegative. There are therefore
no carries among 104-bit digits. Discard digits above the
cutoff, shift the whole integer right by 48, and mask the
low 49 bits of each digit at its original position. This
equals flooring each original digit divided by S: displaced
bits from a following digit begin at bit 104-48=56 and
are erased by the mask. This implements the rounded
probability recurrence using only exact integer operations.

For true anchor depth d within the radius, the lower suffix
coefficients start at C(d). Denote the resulting rounded
lower probability at total child weight z by L_(d,z).
If the child has n_c=1776 observed bits in this program, the mass
of failed child admission at weight z is bounded above by

    binom(n_c,z)*9^(n_c-z)/10^n_c - L_(d,z).

The first term is the exact probability of all child noise
of that weight. The subtracted term is a lower bound on
the admitted portion and the computed difference is
nonnegative.

## Full-change competitor counts

Let M=111, P=996 be the number of parity coordinates, and
Q=2^999 the number of cyclic messages. For e nonzero and x
uniform over the field, the cube difference is

    (x+e)^3+x^3 = e^3 * ((x/e)^2+(x/e)+1).

The linear map z to z^2+z has kernel {0,1}. Its image has
dimension 110 and has trace zero, since Tr(z^2)=Tr(z). It
therefore maps uniformly to that hyperplane. The trace
functional Tr itself is nonzero because 111 is odd and
Tr(1)=1. For a nonzero e the derivative is uniform on an
affine hyperplane with normal the linear functional
Tr(v/e^3). Field multiplication is invertible for every
nonzero multiplier, so the 2^111-1 forms Tr(v/a) for nonzero
a are distinct nonzero linear functionals. Because 3 is
coprime to 2^111-1, e maps bijectively to its nonzero normal.
The linear local parity difference only translates the
hyperplane.

When a parity field is projected to its surviving coordinates,
the affine hyperplane projects uniformly to the whole surviving
space if its normal uses a deleted coordinate. Otherwise the
projection is uniform on the surviving affine hyperplane with
the same nonzero surviving normal. Thus each nonzero character
of a surviving block is realized by at most one nonzero
direction. For a fixed full-change vector e the parity
differences depend on distinct uniform true fields, so their
joint Fourier expansion multiplies the separate uniform or
hyperplane factors.

Fix the true noise and a distance tolerance t. We bound
the number of full-change competitors averaged over the
uniform true message. Expand all parity constraints in binary
characters. Before dividing by the common parity factor 2^P,
the empty parity mask contributes at most the 1995-bit ball
volume V_t.

For a nonempty subset of mask blocks, let k be the number of
full-width blocks in that subset and h=0 or 1 indicate whether
the projected block is selected. These masks are controlled by
(k+h) distinct systematic direction fields, with
c=(k+h)M systematic bits. There are binom(8,k) such subsets for
fixed h. The number of possible mask patterns of weight w
satisfying the nonzero requirement in each selected block is
computed by inclusion-exclusion. For h=0 it is

    sum_(j=0..k) (-1)^(k-j)*binom(k,j)*binom(j*M,w).

For h=1 replace the final binomial by
binom(j*M+108,w)-binom(j*M,w). The projected block's
normal is counted only among its supported masks. For fixed
systematic weight d on the c selected fields there are at most
binom(c,d) choices of their directions. Their masks are distinct,
since all selected nonzero characters determine their direction.
Select the largest binom(c,d) absolute signed parity ball sums
from the mask list with the stated capacities. In the remaining
systematic fields there are at most binom(999-c,d') candidate
patterns of weight d'. Sum over d,d', applying the chosen mask
count at parity radius t-d-d'. Multiplying by binom(8,k) and
summing over k,h gives the upper count.

For a length P parity mask of weight w, the signed sum on the
ball of radius r is the coefficient of X^r in
(1-X)^(w-1)(1+X)^(P-w). This follows by dividing the
generating function for its signed shells by 1-X. Starting
K_(-1)=0, K_0=1, its coefficients satisfy

    (r+1)*K_(r+1)=(P+1-2w)*K_r-(P-r)*K_(r-1),

by differentiating the polynomial. The verifier takes their
absolute values, sorts by magnitude at each radius, and
cumulates counts so ranked partial sums are exact.

All terms carry the density factor 2^-P=8/Q.
Multiply the unscaled sum by eight and call it D_t.
The expected number of full-change competitors is at
most D_t/Q for each fixed true noise. The verifier confirms
that D_t is nondecreasing through the cutoff.

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
resulting first-pair errors. For a full first pair there are eight cell factors, four of width 24, one of width 12 and three of width 1. In the projected first pair the last three factors are 1+X on the unrestricted systematic bit. Each possible
systematic error determines its direction uniquely. The remaining
a-1 systematic and a parity fields have u=(2a-1)*111 bits.
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
this field, here mu=1. The coefficient comes from requiring all k mask
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

If a run starts at field zero, only its first pair is
projected. If a run's changed predecessor set reaches the
anchor parity field, the outside parity coordinates have
three deletions. There are a such starts for length a. For
this case, fill the three absent received symbols with
independent uniform bits for the purpose of the bound.
For every fixed projected noise and candidate, these bits
all agree with the candidate with probability 1/8. On that
event the full and projected distances are equal. A strict
or half-tie indicator in the full space is therefore at
least its projected value times this match indicator.
The unprojected per-noise bound is uniform in the added
bits, so the projected count is at most eight times that
bound. This uses the same D_(a,t) and denominator 2^(p+1),
with p=111*a-3 after projection and u=p+(a-1)*111. No random bits enter the
actual decoder. There are 9-a-1 other starts, disjoint from
the projected fields. These three start counts sum to nine.

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

The checker uses a small fixed set of integer radices a0
for each run length, all stated in RUN_RADICES. At radix a0,
put z=(a0/100)^8. It partitions the first-pair syndrome tuples
by a verified upper bin for H(z). For each component local
polynomial h, the checker evaluates h(z) exactly and rounds it
up to an integer multiple of 2^-76. Let N0 be the product
of the eight rounded numerators for the tuple. For

    k = ceil(log2(N0^8)) - 8*76*8

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
subset of all 222-bit true noise patterns. The checker
verifies the summed bin coefficients are binom(222,d) in a full first pair and binom(219,d) in a projected first pair.

For each bin k the number of true local noises of weight w
is the coefficient of X^w in A_k(X)*(1+X)^u. Representing
polynomials in base 2^(8*stride) packs each coefficient into
one digit. When forming A_k, the stride has room for an
integer of size 2^222; no intermediate coefficient can exceed
the full 222-bit count. When multiplying by the free
binomial factor each output digit through the cutoff is at
most binom(first-pair size+u,w), and the stride is chosen from the largest
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

For lengths 1 through 8 the local weight cutoffs are
66, 92, 119, 145, 172, 199, 226 and 230. For larger local
weights through 230 the checker counts every local pattern
as potentially removed. This is a conservative upper count.
For a given geometry, let S be the total number of affected
positions including the first pair. Convolve the local count
with binom(1995-S,t-w) for the unaffected positions. The
noise probability of each whole pattern at true weight t is
9^(1995-t)/10^1995. Combine the charges for the three
geometry counts. All arithmetic uses the stated upward
rounding and exact remaining operations.

## Occupancy of the actual candidate list

At distance cutoff t<=230 consider all pairs of true word and
received word with true distance at most t. Remove any pair
for which the true word fails the specified anchor or child
admission tests, and any pair for which a single-run competitor
is closer or is tied and lexicographically earlier. The fallback is allowed even when it
would have been removed. Let L_t count removed pairs, Q=2^999,
and V_t be the 1995-bit Hamming ball volume through t. Then
Q*V_t-L_t pairs remain. Two distinct retained candidates at the
same received word must differ in all nine fields, since otherwise
the farther would be removed, or the lexicographically later on a tie.
Their number of ordered distinct pairs at common received words
is at most V_t*D_t by summing over messages and true noises.

Let m_y be the number of retained candidates at a received word
y, and let F_t count received words for which the actual decoder
has some candidate at distance at most t. Every retained candidate
is in the actual list. For every integer k>=1,

    k(k+1)*1[m_y>0] >= 2k*m_y-m_y*(m_y-1),

since the difference for a positive m_y is
(m_y-k)*(m_y-k-1), and both sides vanish at zero. Summing gives

    F_t/Q >= G_t - alpha_t*L_t/Q
    k_t=1+floor(D_t/Q),  alpha_t=2/(k_t+1)
    G_t=floor(V_t*(2*k_t*Q-D_t)/(Q*k_t*(k_t+1))).

Set U_t=max_(i<=t) G_i and choose a maximizing index i(t)<=t.
By monotonicity of F_t we may replace G_t by U_t and charge
alpha_(i(t))*L_(i(t)) in its place. The decoder returns a
candidate of minimum distance in its list including the fallback.
For such a received word at distance t, its contribution to
success is 9^(1995-t)/(Q*10^1995). Thus the success through the
cutoff is obtained by summation by parts with w_t=9^(1995-t)
and w_231=0, by multiplying F_t/Q by (w_t-w_(t+1))/10^1995
and summing. Ignoring success outside the cutoff only decreases
the lower bound. A removed pair with true noise weight v can
be charged only at t with i(t)>=v, and alpha_(i(t))<=alpha_v
because D_t is nondecreasing. Its total charge is at most
alpha_v*w_v/(Q*10^1995) by telescoping.

The error bound therefore consists of the unremoved baseline

    1 - sum_(t<=230) (U_t-U_(t-1))*w_t/10^1995,

plus the proper-run terms for all nine starts, plus the admission
failure charge weighted by alpha at the true total noise weight.
Here U_(-1)=0. Each length-wise charge already sums over all
nine starts and includes the alpha weighting at the true noise
weight. Sum these eight charges once. The overall error upper
bound may count admission and a proper-run removal for one pair.

## Weight-conditioned admission charge

Write d for the true anchor systematic weight on 111 bits.
Write w for its independent parity noise on the other 108 anchor
bits and Z for the child total noise on 1776 bits. The total weight
is d+w+Z, with the three noise blocks mutually independent.
If d>25 the anchor is absent whatever the other noises; its
exact weighted charge is

    sum_(d>25) Pr[Bin(111,1/10)=d]
        * sum_(z<=230-d) Pr[Bin(1884,1/10)=z] * alpha_(d+z).

Use only d<=230 because the weighted charge vanishes otherwise.

For a true anchor depth within the radius, let Z be
the child noise weight above and W0 the anchor parity noise
weight. Write A0 for the distribution of W0, on the 108 surviving
anchor parity coordinates. Condition on d and
z. The child failure event and W0 are independent.
The occupancy weight for this conditional event averages
over its unaffected anchor parity noise:

    G(s) = sum_(w<=230-s) Pr[W0=w] * alpha_(s+w).

For each s the verifier evaluates this rational number
and rounds it upward to scale 10^60. For each d in radius
and each z<=230-d, it uses the exact binomial probability
of d and the upward value for G(d+z), times the upper
child failure mass at z from the exact weight-indexed
suffix calculation. These products are nonnegative.
Summing them and rounding upward to denominator 10^60
gives the child admission charge. The experiment with
part equal to 1 checks the resulting integer in joint.json
by deriving all of its terms. The experiment with part
equal to 3 uses that verified upper charge together with
the exact anchor and checked proper-run contributions.

The exact probability of successful repetition decoding is
sum_(j=0..2) binom(5,j)*9^(5-j)/10^5=99144/100000.
The two codes use disjoint independent noise coordinates and
independent uniform information. If the certified cyclic error
upper bound is E, the total error is therefore at most
1-(1-E)*99144/100000.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The eighteen received fields are copied to padded little-endian
two-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
The working arrays are consecutive from word 48.
The fixed data include reduction and square tables, the Golay and repetition mapping and candidate tables, and four additional
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
- QR_ROW_BASE, one sorted 2^23-entry list per representative,
  with each low 48-bit raw noise word and its cost above it.
- BUDGET_BASE, one word per anchor depth 0 through 25 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (widths 12 and 1) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*(3301)+r for every type p and state r.

Construction.py, qr.py and dynamics.py specify every fixed datum. The exclusive end of RAM is the value MEMORY_END in
construction.py, within the model's 2^32-word limit. All table addresses are distinct and the
program's temporary words precede these tables.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All 64 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use two accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=2 words and h=47 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^111. This fits in one word since the
largest possible degree is 14. To multiply a field
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
c*X^111 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=14 byte slots, emit a
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
Move R registers to zero. For each of n0=28 nibble
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

Extract a literal piece of width at most 64 by loading its
first word and shifting by its offset. If it crosses a word,
load the next word, shift by the complementary amount and OR.
Mask the result to the literal width. At most six primitives
suffice. Emit this operation and a store for each of 36 padded
Y and P words using the encoder's offsets, padding P0 after
its 108 surviving bits. Form the all-zero incumbent distance
by loading, popcounting and summing those 36 received words
and initialize all eighteen B words to zero. Copy Y0 to X0
in four primitives. Set Depth to zero and store NEXT[0]=0.
No initialization or output is inside the anchor traversal.

For packing a field from known short pieces, clear two
registers by free moves and insert a piece with a left shift
and OR. At a crossing shift its original value right and OR
into the next word. Store the two accumulators. All indices,
widths and crossings are literal.

At output gather the nine incumbent systematic fields into
the sixteen output words. Add the bit at input position 999 and
the four bits at positions 1996 through 1999 and put the
comparison sum>=3 at output position 999. The 24 high output
bits are zero. Each output word uses at most four literal
pieces and a store, with at most eight primitives per piece.
The five extra bit extractions, summation and comparison use
at most 20 primitives before inserting the Boolean. Charging
252 for the received field copies, 108 for the zero distance,
18 for clearing the incumbent, 5 for anchor initialization,
20 for majority and 16*33 for output gives at most 931
nontraversal primitives.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=111, and Depth is a persistent register.

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
        branch if done to output packing
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
N=sum_(i<=R) binom(111,i) visited anchors. These events
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
This uses one cube and at most 250 further primitives:
the eight cells use at most 8*(6+2+4), and finishing AP
uses 6, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 64 cells in field and within-field order.
Its frame word slots at offsets 0..5 are Base, U, Rem,
Dist, Rank and Perm. The literal SELECTED slot holds
this cell's selected systematic value. Parent continuation
means the previous block's CELL_NEXT, or ANCHOR_EXTEND
for the first block.

For width m=12 or 1 use the following emitted program.

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
        j := mul(p,3301)
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
        needed := add(cc,64-DepthCell)
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
        run FIELD after each completed nonfinal child field
        run FINAL if this is the literal final cell
        jump to the next CELL_INIT or, at the final cell,
            back to its CELL_NEXT

Each nested load or arithmetic expression is a separate primitive.
Initialization uses at most 36. The unique stopping attempt
at a visited frame costs at most 11. Each successful attempt
costs at most 28 including its jump, plus FIELD or FINAL where applicable.

At a width-24 cell change the initialization to use the same three
extracts, the QR_PARITY_BASE lookup and the same two XORs to form
s. Then emit

    addr := add(QR_ORBIT_BASE,s)
    label := load(addr)
    orbit := and(label,511)
    perm := shr(label,9)
    row := mul(orbit,2^23)
    base := add(QR_ROW_BASE,row)
    pos := mul(perm,6*256)
    pb := add(QR_PERM_BASE,pos)
    store(Base,base)
    store(Perm,pb)
    store(U,u)
    idx := add(PROFILE_BASE[24],orbit)
    p := load(idx)
    j := mul(p,3301)
    j := add(j,Remaining)
    idx := add(CAP_OFFSETS[DepthCell-1],j)
    effective := load(idx)
    store(Rem,effective)
    store(Dist,Distance)
    store(Rank,0)
    fall through to CELL_NEXT

The two loop tests and ordering in CELL_NEXT are the same, using
rank limit 2^23 and cc=shr(record,48). After the successful
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
The QR initialization uses at most 43 primitives, including
its three initial extracts, and its single failed attempt at
most 11. Charge at most 67 for these together and at most
80 for each successful attempt including its jump. Each field boundary falls at an ordinary width-1 cell. All stages use the
six frame words and fewer than 20 scratch registers alongside
the persistent registers specified above.

When an accepted child returns, the parent still has its
Base, U, saved budget and distance and next rank in RAM.
The frame stops when it runs out of entries or finds the
first over-budget record. The costs are sorted, so further
entries would also fail. Each full tuple in the specified
list is therefore visited exactly once.

At FIELD and FINAL, pack the eight selected cell slices from
the current completed field into its two padded Xj words,
with at most 8*(1+4)+2=42 primitives. Compute Cj=Xj^3.
At FIELD this cube supplies the next field; its total routine
charge is at most CubeOps+180. At FINAL add the distance of
the anchor parity to Distance by loading AP and C8 round by
round over two words, XORing, popcounting and adding.
Before the second popcount, AND the XOR with 2^44-1, retaining
only the 44 transmitted bits of that anchor word. This adds one
primitive within the 210 noncube allowance.
The retained bits of AP+C8 equal the transmitted-coordinate
projection of P0+L(X0)+C8. Distance already includes
the anchor systematic distance and the eight visited child pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all eighteen systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 210 primitives.
## Exact worst-case operation bound

Number the 64 literal cells by i=1..64. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 67 at width 24 and 47 otherwise, and g_i be 80 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
T_i is CubeOps+180 at i=8,16,24,32,40,48,56 and CubeOps+210 at i=64.
It is zero at all other indices.

For all possible states r=0..3300 define in decreasing cell order

    F_65(r)=0
    U_(i,p)(s)=sum_(c<=s-(64-i))
                 A_(p,c) * (g_i + T_i + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,3300]: U_(i,p)(s)<=F_i(r)-h_i}.

The nonnegative recurrence is nondecreasing in r by induction
and U is nondecreasing in s. The maximum defining E exists
because r itself is allowed. Its search returns the rightmost
valid index, including equality. These are the fixed CAP table entries.
For any actual state r and profile p, the accepted attempts
are exactly the rows up through their recorded cost threshold
E_(i,p)(r)-(64-i). Their continuation states E_(i,p)(r)-c
belong to [0,3300] because costs are positive and successful
attempts satisfy the reserve. At a successful attempt the
remaining program cost is bounded by g_i+T_i+F_(i+1) at that state.
Consequently the actual entire cell traversal costs at most
h_i+U_(i,p)(E_(i,p)(r))<=F_i(r). This proves the bound by
downward induction without distributional hypotheses.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 250+CubeOps. Hence a worst-case
all-input integer bound for the single cyclic traversal and its layout is

    W = 1000 + sum_(d=0..25) binom(111,d) *
                 (44 + 250 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.

## Exact score computations

The experiments verify the QR spectra, field checks, orbit
count and worst-case primitive bound in part 0. Part 1
certifies the integer upper child charge in joint.json at
denominator 10^60. Parts 2 and 3 compute the length-wise
proper-run ceilings, rounded upward to denominator 10^18,
for lengths 1 through 5 and 6 through 8, respectively.
Together they verify all eight entries in losses.json. In
part 3, combine these verified ceilings with the exact
baseline and anchor terms and the verified child charge
before applying the independent majority success probability.

    Unremoved cyclic baseline error <= 61542415046278112/10^18
    Proper-run charges for all starts by length <= (96721937552531, 2229493958972482, 5823632443536058, 8266993049988186, 5939192471722964, 3182580798770005, 856713942379716, 27832504256259)/10^18
    Anchor admission charge <= 20739251387552/10^18
    Child admission charge <= 27446739619985555/10^18
    Cyclic error <= 115433055024829419/10^18
    B = 123004948073816879/1000000000000000000
    W = 123188796072353641708753665975402158045817211354180
    W^500 <= 2^83199
    Exclusive RAM end = 3705364644.

The score rounds the final exact product upper bound upward.
The component display rounds separately.
