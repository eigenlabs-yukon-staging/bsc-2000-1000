# Shortened cyclic fields with exact weight-indexed admission

## Encoder and restricted decoder

Work in F2[X]/(P) with P=X^334+X^5+X^2+X+1 and the polynomial
basis, low coefficient first. The prime divisors of 334 are 2 and 167. The experiment checks
X^(2^334)=X modulo P and the two coprimality tests with
X^(2^2)+X and X^(2^167)+X required by Rabin. Let U consist of the field
values whose two highest coefficients are zero.

Split the 1000 message bits into a 332-bit field value b0 in U
and two 334-bit field values b1 and b2. Append three cyclic parity fields given by

    p_i = L(b_i) + (b_(i-1))^3.

Parity field p1 has its two highest bits omitted;
all others transmit all 334 bits. Thus there are 1000 systematic
and 1000 parity channel bits. The encoder is injective by its
systematic prefix. Construction.py fixes all layouts.

Each L consists of thirteen width-24, one width-12, two
width-4, and two width-1 identity cells in increasing bit order.
For a width-24 input u use low 23 parity bits equal to the
remainder of u*X^23 modulo the degree-23 binary G=9205487.
Set its high parity bit to the parity of the weight of u and
its remainder. The width-12 Golay cell uses the same rule
with degree 11 and polynomial 1+X^2+X^4+X^5+X^6+X^10+X^11.
The width-4 cell has columns 14,13,11,7 and the last two
width-1 parity cells are identities. All maps are linear.
In particular the high two coordinates of L(b0) are zero.

The decoder visits each anchor x0 in U within 71 systematic
errors of Y0. At depth d (the anchor systematic distance), use
C(d) from entry d of construction.py BUDGETS. For a fixed anchor, visit child fields 1 and 2
in order, and their literal cells in the specified order. A paired
cell of width m has systematic observation u and parity observation
v adjusted by XOR with the already chosen predecessor cube.
Write s=v+L_m(u) for the syndrome. For a possible systematic
error e its paired distance is

    D_s(e) = wt(e) + wt(L_m(e)+s)
    A_s(e) = 9^(2m-D_s(e))
    S_s = sum_e A_s(e)
    c_s(e) = ceil(16 log2(S_s/A_s(e))).

For the two unpaired width-1 cells of field 1, the observation
is only its systematic bit. Use D(e)=e, A(e)=9^(1-e) and
S=10 in the same cost formula. All costs are positive integers,
since all rows have several positive entries. The preloaded cost
is exactly the bit length of (S^16-1)//A^16, since
S^16<=2^h*A^16 if and only if that quotient is below 2^h.
Width-24 candidate rows use the first 2^23 of the 2^24
possible raw noise patterns in a syndrome orbit representative.
The fixed permutation described below maps them back to the
actual syndrome. All candidates are retained at other widths.
The rows are sorted by distance and then by their stated
numeric error order, so cost is nondecreasing.

Depth-first search the 36 child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(36-i), where r is
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
and the anchor. Cell i has b_i observed bits, with b_i=2m
for a paired width-m cell and b_i=1 for a systematic-only
cell. Let H_(i,p,c,w) be the number of retained true noise
patterns of weight w and nominal cost c across all
syndromes of its profile type p. Multiply the retained
count for one syndrome by its profile multiplicity.
The QR orbit permutation preserves total weight and cost,
even if truncation cuts the last retained weight tie.
In a systematic-only cell there is one type, and one
pattern each at weights zero and one with their stated
costs. The channel mass for a term (p,c,w) is

    m_(i,p,c,w) = H_(i,p,c,w)*9^(b_i-w)/10^b_i.

For each i, use the same effective cap E_(i,p)(r) and
remaining state r as the decoder. With s=E_(i,p)(r),
a true pattern of cost c passes precisely when it is retained in
its row and c<=s-(N-i), where N is the fixed number of
child cells in the program. The next state is s-c.
For the true path this follows even in a shortened cell,
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
through 238 are required, because the scoring contribution
in this proof only uses total noise weights up to 238.

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
noise patterns at each total weight v through 238.
Part 2 computes and checks all the resulting counts in misses.json.

## Competitors changing all three fields

For a full degree-334 field, write q=2^334. For e nonzero,
the cube difference is

    (x+e)^3+x^3 = e^3 + e*x^2 + e^2*x.

The map z -> z^2+z has kernel {0,1} and image the kernel
of the absolute trace. The trace is an F2-linear nonzero
functional by the polynomial root bound. Thus for x uniform
over the full field, the parity difference is uniform on an
affine hyperplane. Its normal has entries Tr(e^(-3)*X^j).
The pairing (h,z) -> Tr(hz) is nondegenerate because a nonzero
h permutes the field and trace is nonzero. Thus a normal
determines the field multiplier e^(-3). The cubic root bound gives at most three directions per normal.

For parity p1 the predecessor x lies in U, and only the low
332 coordinates of its difference are transmitted. Any mask
on the transmitted coordinates lifts to a unique multiplier h
for the full trace pairing, by extension as a functional
on the 334 basis coordinates. For a nonzero mask h is nonzero.
It is constant on the derivative image of U exactly when

    L_e^*(h) = (h*e)^(2^333) + h*e^2

belongs to U-perp, the orthogonal complement under trace.
Indeed the displayed adjoint follows because
Tr(h*e*x^2)=Tr((h*e)^(2^333)*x). The space U-perp has
four elements. Fix any z from it and square the equation
L_e^*(h)=z to obtain

    h*e + h^2*e^4 = z^2.

This is a nonzero degree-four polynomial in e, with at most
four roots. For each nonzero mask there are therefore at most
16 possible nonzero directions for which it occurs. No
assumption about the actual projection rank is needed for this
per-mask bound.

When all fields change, fix the differences e_i. Parity
differences in the three transmitted fields are independent because
they use distinct true predecessor fields, which are independent
under the uniform systematic prior. Within each field the
distribution is uniform on a translate of a binary linear image.
If that field has parity length n, its probability at a vector
is 2^-n times a sum of characters over all annihilating masks,
with translation signs. For example on an affine support t+V
the sum of the characters indexed by V-perp is |V-perp| at x
in the support and zero outside; dividing by 2^n gives 1/|V|.
The linear term L(e_i) is just a translation, including after
puncturing. Multiplying all three expansions gives a factor 1/Q
with Q=2^1000, and every term is an empty contribution or
a nonempty combined parity mask.

For a mask of weight w>0 in the 1000 parity bits, the absolute
signed radius-l ball sum is abs(K_l(w)), where K_l(w) is the
coefficient of X^l in
(1-X)^(w-1)*(1+X)^(1000-w). To see this, the signed
per-shell polynomial has exponent w on (1-X), and summing
prefixes divides by 1-X. Let g=1001-2w.
Differentiating gives the exact coefficient recurrence
(l+1)K_(l+1) = g*K_l - (1000-l)*K_(l-1),
with K_(-1)=0 and K_0=1. The checker verifies divisibility.

Suppose the term selects k of the two ordinary parity masks.
Optionally it also selects a nonzero mask on p1, indicated by
h=1 (otherwise h=0). There are binom(2,k) choices of ordinary
positions. A combined weight-w mask appears in at most the
following counts across the corresponding predecessor directions:

    h=0: 3^k*[Z^w] ((1+Z)^334-1)^k
    h=1: 16*3^k*[Z^w] ((1+Z)^334-1)^k * ((1+Z)^332-1).

These follow from ignoring zero masks at selected positions
and the three or sixteen direction preimages. The calculation
evaluates those coefficients by inclusion and exclusion.
There are n=334k+332h systematic bits at the selected
predecessor positions, with the other 1000-n systematic bits
unrestricted by these masks. For fixed true noise and
combined candidate systematic error weight a on selected
predecessors, there are at most binom(n,a) choices of their
directions, since each possible systematic error determines them.
For a fixed nonzero direction e in U the linear derivative on U
has kernel {0,e} and image dimension 331. Projection can
remove at most two dimensions, so its rank is at least 329
on the 332 parity coordinates. There are consequently at most
7 nonzero annihilating masks for this direction. Ordinary selected
predecessors each have only one nonzero mask. Thus each tuple
of directions contributes at most one combined mask when h=0
or at most seven when h=1.

For each l rank abs(K_l(w)) by the listed per-mask multiplicities.
Take the sum of the largest binom(n,a) values when h=0, or
the largest 7*binom(n,a) values when h=1. This accounts for
all (direction tuple, mask) pairs at that systematic weight.
The checker fills groups of the largest magnitudes including
a partial last group where necessary.

For candidate systematic error weight b on the other predecessors
multiply by binom(1000-n,b). Set l=t-a-b and sum over
a,b with l>=0, then all nonempty mask selections. The empty
selection contributes the full 2000-bit Hamming ball volume at
radius t by Vandermonde convolution over the 1000 systematic and
1000 parity bits. Call the resulting integer D_fixed(t). The expected
number of full-change competitors at distance at most t, for every
fixed actual noise, is bounded by D_fixed(t)/Q. The checker verifies
the bound's nondecreasing values through t=238.

## Averaged full-change ball intersections

The occupancy calculation uses the total number of ordered
full-change pairs with a common received word within distance
t of both codewords. For cutoffs 198 through 238 we bound this
count directly by summing ball intersections over all actual
true noise patterns in the ball.

Let M=334, G=3, K=1000 and Q=2^K. Write J_t(d) for the
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

Let A=332 be the shortened dimension and let
delta=M-A=2. Write N_(k,s)(d) for the number of
weight-d patterns with k ordinary nonzero blocks and s
nonzero shortened blocks, where s is either zero or one.
Inclusion-exclusion on the ordinary blocks gives

    N_(k,s)(d)=sum_(j=0..k) (-1)^(k-j)*binom(k,j)*F_s(j,d),
    F_0(j,d)=binom(j*M,d),
    F_1(j,d)=binom(j*M+A,d)-binom(j*M,d).

The subtraction in F_1 requires the shortened block to
be nonzero. The single shortened systematic field and
single shortened parity field belong to the same affine
difference factor. Fix a subset selecting k ordinary
factors and s shortened factors for nonzero parity masks.
The mask capacity at weight w is at most

    mu^k*(4*2^delta)^s*N_(k,s)(w),

where mu is the ordinary cubic-normal multiplicity bound.
The shortened mask capacity accounts for every lift of a
projected parity normal and the degree-four preimage bound
proved above. For any one shortened direction the
annihilator contains at most eta=2^(delta+1)-1 nonzero
masks. Keep capacities by exact mask weight. At complement
systematic weight b, there are N_(G-1-k,1-s)(b) patterns,
and selected fields supply at most eta^s*N_(k,s)(d_s-b)
terms. Sum that many largest upper terms within the
exact per-weight capacities, multiply by the complement count, and
sum over b. Multiply by binom(G-1,k) for the possible
ordinary subsets. Sum over all nonempty selected subsets
and d_s=G..2t to obtain the integer correction R_t.

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
Parts 1 and 2 together check every declared integer in pairs.json.

## Proper-run competitors

A proper nonempty changed set separates uniquely into cyclic
runs with disjoint affected coordinates. A length-a run starting
at i changes its a systematic fields and the parity fields i
through i+a. A proper competitor with distance no greater than
the true one forces a single run with distance no greater,
since the affected differences add. Charge all ties here.

The first pair (systematic field i and parity field i) changes by
direction e and transmitted L(e), as its predecessor is unchanged.
Each paired syndrome is invariant. For fixed true local noise let
H_s(Z)=sum_d h_s(d) Z^d describe the possible errors at a paired
cell. If i=0, its two transmitted parity bits with fixed zero
systematic coordinates also stay unchanged. For each their fixed
actual noise bit s yields the enumerator Z^s. If i=1, its two
unpaired systematic cells have unrestricted error patterns with
enumerator 1+Z and only one possible syndrome class. All cases
therefore have eighteen small factors in the first-pair polynomial H.

Let p be the total length of parity fields i+1 through i+a,
let u be p plus the systematic length of the other changed fields,
and let n be u plus both lengths in the first pair. The parity
differences among these p bits are independent for fixed directions.
If no changed predecessor is b0, their point mass is at most
2^a/2^p. If b0 changes, its direction e is a nonzero element
of U. The map e*x^2+e^2*x has kernel {0,e} on U, so
its image on the 332-dimensional domain has dimension 331.
Puncturing can lose at most two independent coordinates and
leaves image dimension at least 329. Its point mass on the
332 transmitted coordinates is at most 8/2^332. Thus in either
case the product mass on the p bits is bounded by

    2^(a+2*j-p),  j=1 exactly when b0 lies in the run.

For fixed exact true noise, the conditional expected number of
single-run competitors of local distance at most w, the actual
weight on the n affected bits, is bounded by this mass times
the coefficient prefix through w of (1+Z)^u*H(Z). A potential
systematic error fixes all run directions, and the invariant
syndrome enumerators count their possible first-pair errors.
The point mass is used only on actual nonzero directions
before enlarging the resulting coefficient count to all patterns.

For any rational r from {28/100,30/100,...,54/100}, evaluate at
Z=r^2. The prefix is at most r^(-2w)*(1+r^2)^u*H(r^2).
The conditional probability of any such competitor is at most
min(1,x) for the expected-count bound x, and min(1,x)<=sqrt(x).
Consequently this probability is at most

    2^(-floor((p-a-2*j)/2)) * r^(-w)
      * (1+r^2)^ceil(u/2) * product_cells sqrt(H_s(r^2)).

Both rounding directions increase the bound. For an evaluated
cell enumerator H_s(r^2)=A/B, the checker rounds its root up
to a multiple of 2^-64 by computing j0=isqrt(floor(2^128 A/B))
and adding one unless j0^2 B=2^128 A. For each first-pair
cell multiply this root numerator by h_s(d) and sum over s
at each true weight d. This covers exactly all channel noise
patterns, including the fixed-parity and unpaired cases above.

Form the product of the eighteen polynomials so obtained and
(1+Z)^u. At each local weight w multiply by the rational
common factors above, dividing by the eighteen powers of 2^64
for the roots. This bounds the sum of the individual probability
bounds over true local noise patterns of that weight. Minimize
the result at each w over the stated grid and the trivial
binom(n,w) bound. Each candidate is uniform for every fixed
outside noise of that weight class, so the minimum before outside
expansion is valid. The local sum is then multiplied by
binom(2000-n,v-w).
Sum over w to get a bound J_(i,a,v) on each average removal count
for start i, length a and total true weight v. The checker
groups exactly identical geometries by their multiplicities.

Each syndrome enumerator has at most 2^24 patterns. Each
pattern contributes at most one at r^2<=1, so H_s(r^2)
is at most 2^24. Its square root scaled by 2^64 is at most
the integer 2^76, and its upward ceiling is also at most
2^76. Each product of paired components sums to at most
2^(668+18*76); the free factor
adds at most 2^1002. Thus 4000-bit digits suffice throughout
the packed polynomial products. All further arithmetic is rational.

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

At cutoff t<=238 let V_t be the volume of the 2000-bit ball,
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

    sum_(t=0..238) (G_t-G_(t-1))*9^(2000-t)/10^2000.

Part 3 checks the claimed error bound from this exact lower
success bound and the validated witness counts.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The six received fields are copied to padded little-endian
six-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed data include reduction and square tables, the Golay, Hamming and repetition mapping and candidate tables, the two-entry systematic-only row and four additional
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
- BUDGET_BASE, one word per anchor depth 0 through 71 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (widths 12, 4 and 1) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*4601+r for every type p and state r.
  The single systematic-only profile type p=0 has no type table.

Construction.py, qr.py and dynamics.py specify the original RAM
regions, including the fixed admission tables. Accelerator.py adds
multiplication-row scratch words starting at their exclusive end,
MEMORY_END. The fixed wide-reduction and wide-square tables follow
that scratch range. Their exclusive end FAST_END fits the model's
2^32-word limit, as checked in the declared experiment. Fixed and
working regions are disjoint.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All 36 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube macro uses six vectors of 6 words each and at most six
additional scratch registers. Together with the four persistent
registers this uses at most 46. Other routines' scratch
registers are reused where their lifetimes do not overlap.

### Exact field cube

Write k=6 for the field word count and h=14 for the
number of meaningful high-word bits. Let P be the field modulus
and P_low=P+X^334. The following three stages are the emitted
cube_schedule in accelerator.py. It replaces only the cube macro
in the rest of the stated decoder. Every cube yields the same
field value, so the transmitted code, candidate admission, costs,
choice of candidate, and output bits in the error proof are
unchanged. The fixed budget and cap tables remain the ones
specified by dynamics.py.

First read the k padded source words into registers. Partition
the source coefficients into successive 22-bit chunks, with a
short final chunk if needed. For each chunk index j and possible
value d, the wide-square table stores k words of (d X^(22*j))^2
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
at most 16. One address add and one load retrieve the word
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
primitive instructions gives CubeOps=3405 for all sources.
Free register moves are excluded under the model. The experiment
checks the exact opcode counts and replays eight representative
field values on the unchanged organizer counter. The replay
shows consistency of the fixed schedule and field arithmetic.
The identity of the cube for every field value follows from
the linear square, modular shift and Horner arguments above.

### Bit pieces and initial state

To extract a literal bit piece of width 1 through 64, load its
starting word, and right shift by the within-word offset if
nonzero. For a crossing, also load the next word, shift left
by the complementary offset, and OR. If the part is shorter
than 64 bits, mask it. This costs at most six operations.
Use this for the 36 padded words of Y and P, using the
transmitted widths in construction.py, then store each.
The shortened Y0 and P1 have two extra zero padding bits.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 108. Store zero in the 18
padded incumbent words. Copy Y0 to X0 in twelve primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear six
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the six
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. For output use field widths 332,334,334. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=332, and Depth is a persistent register.

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
N=sum_(i<=R) binom(332,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base (QR_PARITY_BASE at width 24),
load the parity image, and pack
these images into a six-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining by adding Depth to BUDGET_BASE and loading
that word. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 250 further primitives:
the 18 cells use at most 18*(6+2+4), and finishing AP
uses 18, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 36 cells in field and within-field order.
Its frame word slots at offsets 0..5 are Base, U, Rem,
Dist, Rank and Perm. The literal SELECTED slot holds
this cell's selected systematic value. Parent continuation
means the previous block's CELL_NEXT, or ANCHOR_EXTEND
for the first block.

For a paired width m=12,4,1 use the following emitted program.

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
        j := mul(p,4601)
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
        needed := add(cc,36-DepthCell)
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
        run FIELD after completing the nonfinal field
        run FINAL if this is the literal final cell
        jump to the next CELL_INIT or, at the final cell,
            back to its CELL_NEXT

Each nested load or arithmetic expression is a separate primitive.
Initialization uses at most 36. The unique stopping attempt
at a visited frame costs at most 11. Each successful attempt
costs at most 28 including its jump, plus FIELD or FINAL where applicable.

For the two systematic-only cells of field 1, use width one.
Instead of the ordinary initialization, extract the systematic
bit and perform the same five stores, using the literal
SYS_ONLY_BASE as Base. The profile type is always zero.
Obtain effective with an add of Remaining to the literal
CAP_OFFSETS[DepthCell-1], then load; store it in Rem. The row has exactly two candidates,
with the indicated one-bit distance and normalized cost.
CELL_NEXT and its width-one rank stop are unchanged.

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
    j := mul(p,4601)
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
80 for each successful attempt including its jump. Each field boundary falls at an width-1 cell. All stages use the
six frame words and fewer than 20 scratch registers alongside
the persistent registers specified above.

When an accepted child returns, the parent still has its
Base, U, saved budget and distance and next rank in RAM.
The frame stops when it runs out of entries or finds the
first over-budget record. The costs are sorted, so further
entries would also fail. Each full tuple in the specified
list is therefore visited exactly once.

At FIELD and FINAL, pack the eighteen selected cell slices from
the current completed field into its six padded Xj words,
with at most 18*(1+4)+6=96 primitives. Compute Cj=Xj^3.
At FIELD this cube supplies the next field; its total routine
charge is at most CubeOps+180. At FINAL add the distance of
the anchor parity to Distance by loading AP and C2 round by
round over six words, XORing, popcounting and adding.
The residual AP+C2 is P0+L(X0)+C2. Distance already includes
the anchor systematic distance and both visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all eighteen systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 200 primitives.
## Exact worst-case operation bound

Number the 36 literal cells by i=1..36. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 67 at width 24 and 47 otherwise, and g_i be 80 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
To reproduce the fixed admission caps, use the original cube
cost parameter CubeOld=5037 in the cap definition.
Let T_i^0 be CubeOld+180 at i=18, CubeOld+200 at i=36,
and zero elsewhere. Define T_i^* the same way with CubeOps
in place of CubeOld.

The unchanged cap recurrence has states r=0..4600.

    F_{37}(r)=0
    U_(i,p)(s)=sum_(c<=s-(36-i))
                  A_(p,c) * (g_i + T_i^0 + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,4600]: U_(i,p)(s)<=F_i(r)-h_i}.

All terms are nonnegative, so F and U are nondecreasing
by downward induction. Each maximum exists since s=r is
eligible. Dynamics.py derives these exact cap values, using
the original cube cost solely for this data definition.
This is also the cap table used by the true-path probability
calculation. No cap is changed by accelerating the cube.

For the actual operation bound define G_{37}(r)=0.
For each cell i, set

    U^*_(i,p)(s)=sum_(c<=s-(36-i))
                   A_(p,c) * (g_i + T_i^* + G_(i+1)(s-c))
    G_i(r)=h_i+max_p U^*_(i,p)(E_(i,p)(r)).

For any actual state r and profile p, the emitted list traversal
accepts exactly the entries through its fixed cap minus the
reserve. The resulting states E_(i,p)(r)-c stay in [0,4600].
By downward induction the accepted attempts cost at most the
terms in U^*, and h_i covers the setup and the stopping attempt.
Thus G bounds every execution. This direct calculation uses
the old admission caps and the new primitive work, without a
worst-case claim based on inputs sampled from the channel.

The anchor subset traversal uses at most 44 primitives per visit,
and its parity residual setup and cube at most 250+CubeOps.
Together with fixed initialization and output the resulting
all-input bound is

    W = 1000 + sum_(d=0..71) binom(332,d) *
                 (44 + 250 + CubeOps + G_1(C(d))).

The declared experiment regenerates the fixed caps, computes this
recurrence, and checks the exact exponent inequality.

## Scores

The exact counted primitive bound is W=2223356814838923194034268079322755995552747085883640128816673200260242659635641535567724982759906, with
W^11861 <= 2^3796206. The verifier checks this
inequality and the exclusive RAM address limit.

The exact success calculation above gives error at most
B=27587659089315809896980486855893974712111/500000000000000000000000000000000000000000. The declared experiment checks that rational inequality
using integer or fraction arithmetic throughout.
