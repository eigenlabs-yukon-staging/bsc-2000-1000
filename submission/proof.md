# Four independent five-field components with exact weight-indexed admission

## Encoder and decoder

Use the polynomial basis of F2[X]/(P), where
P=X^50+X^9+X^2+X+1. The declared field check applies
the Rabin criterion with prime divisors 2 and 5 of 50.
Split the 1000-bit message into four 250-bit parts, each split
into five 50-bit words b0 through b4 in increasing bit order.
Within each part define p_i=L(b_i)+b_(i-1)^3 with cyclic
indices modulo five. Transmit the four parts' systematic words
in order, then their parity words in the same part and field order.
Each part is a 250-bit systematic encoder with 500 channel
coordinates. The 2000-bit encoder is injective by its systematic
prefix. The decoder below processes each component on its own
500 received coordinates and concatenates their four 250-bit
results. Reused scratch and tables do not introduce dependence
between decoded component values.

In each field L is block diagonal with two 24-bit cells and two
width-1 repetition cells in increasing bit order. For a width-24
input u, the low 23 parity bits are the binary polynomial
remainder of u*X^23 modulo the degree-23 G=9205487 (in
low-coefficient-first order). The high parity bit is the parity
of the sum of Hamming weights of u and the low parity bits.
The width-1 parity map is the identity. These maps are linear. Construction.py fixes the integer layout.

The decoder enumerates all anchors x0 of distance d<=11
from received Y0. For each use entry C(d) at index d of construction.py BUDGETS as the initial
remaining-cost state. Visit the child fields j=1 through 4 in order. In each field,
form the local systematic and adjusted parity observations u,v
from received Yj and Pj+x_(j-1)^3. At width m define the local
syndrome s=v+L_m(u). For each error e on the systematic
part, the paired distance and unnormalized weight are

    d_s(e)=wt(e)+wt(L_m(e)+s),
    A_s(e)=9^(2m-d_s(e)),   S_s=sum_e A_s(e).

Assign cost c_s(e)=ceil(16 log2(S_s/A_s(e))). The preloaded
integer cost is bit_length((S_s^16-1)//A_s(e)^16), since
the quotient is below 2^c exactly when S_s^16<=2^c*A_s(e)^16.
Costs are positive because every complete row has more than
one positive term.

All candidates of width 1 appear in fixed rows
ordered by (distance, numerical systematic error). In width 24
use a row of the first 2^23 noise patterns, through an explicit
syndrome orbit map and coordinate permutation defined below.
This preserves the distances and the syndrome. The selected
cell is u+e, where e is the systematic part of its noise pattern
in original coordinates. All rows are in nondecreasing distance
and thus nondecreasing cost.

Depth-first search the 16 child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(16-i), where r is
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
scale S=2^48. Round each local mass down to
floor(S*m_(i,p,c,w)). The empty suffix polynomial has
coefficient S at degree zero. At a step multiply a lower
suffix coefficient by the appropriate local mass integer
and sum contributions, then floor the resulting coefficient
divided by S. Nonnegativity gives, by induction, lower
bounds on S times all the exact P coefficients. Only degrees
through 65 are required, because the scoring contribution
in this proof only uses total noise weights up to 65.

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

For a true anchor depth d within the radius, start the suffix
at C(d). Let l_(d,z) be its lower integer probability numerator
at child weight z, with denominator S=2^48.
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
noise patterns at each total weight v through 65.
Part 2 computes and checks all the resulting counts in misses.json.

## Competitors changing all five fields

Write q=2^50 and Q=q^5. For nonzero e, the cube
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
annihilate trace on a spanning set. Each nonzero normal mask
determines at most three nonzero directions e, because after
inverting the trace coordinate map and the field inverse,
the remaining equation is cubic over a field.

A full-change competitor is specified by nonzero
directions e_i in all five fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the five translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 250
parity bits can be evaluated by expanding the
five affine indicators. Multiplication by their
common point mass 32/Q leaves 1/Q times a sum
of 32 terms. The empty mask term is the
parity ball volume. Each nonempty term has an
absolute bound equal to the absolute signed sum
of the parity ball against the combined mask of
its selected parity components.

For a 250-bit mask of weight w>0 the signed ball
sum is, up to translation sign, the coefficient
K_l(w) of X^l in (1-X)^(w-1)*(1+X)^(250-w).
Indeed the polynomial of shell signs has the same
form with exponent w on (1-X), and summing the
prefix divides by 1-X. Let N=249 and
g=N-2*(w-1). Differentiating this polynomial F gives
(1-X^2)F'=(g-NX)F. The coefficient recurrence is

    (l+1)*K_(l+1)=g*K_l-(N-l+1)*K_(l-1),

with K_0=1 and K_(-1)=0. The calculation checks each
division exactly.

In a term that selects k masks, their combined
weight w has at most the following multiplicity
across ordered directions on the contributing
systematic fields,

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(50*h,w).

The sum before multiplying by 3^k counts weight-w combinations
with no zero mask on any of the k blocks. This is
inclusion and exclusion from all binary weight-w
words on the k blocks. Each mask has at most three
direction preimages, giving the factor 3^k. For a given l and k define
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
most binom(50*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(50*(5-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(50*(5-k),b) *
    H_k(t-a-b, binom(50*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(5,k) selections of masks, and add
the 500-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_fixed(t). The expected
number of full-change competitors within total
t for each fixed true noise is at most D_fixed(t)/Q.
The checker verifies D_fixed(t) is nondecreasing through t=65.

## Averaged full-change ball intersections

The occupancy calculation uses the total number of ordered
full-change pairs with a common received word within distance
t of both codewords. For cutoffs 1 through 65 we bound this
count directly by summing ball intersections over all actual
true noise patterns in the ball.

Let M=50, G=5, K=250 and Q=2^K. Write J_t(d) for the
intersection size of two radius-t balls at distance d.
On the d differing coordinates choose a differences from
the first word and d-a from the second. The remaining equal
coordinates can contribute at most t-max(a,d-a) differences.
Therefore

    J_t(d)=sum_(a=0..d) binom(d,a) *
                 sum_(b=0..t-max(a,d-a)) binom(500-d,b).

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

Write N_k(d) for the number of d-weight patterns on k
M-bit blocks with each block nonzero. Inclusion-exclusion gives

    N_k(d)=sum_(i=0..k) (-1)^(k-i)*binom(k,i)*binom(i*M,d).

In particular N_0(0)=1 and N_0(d)=0 for d>0. Fix a subset
of k selected parity blocks. Its masks have distinct
predecessor systematic fields. The total weight-w mask
capacity among actual direction tuples is at most
mu^k*N_k(w), where mu=3 in this even degree, by the cubic-normal
bound. Keep the capacities at each exact mask weight. For complement
systematic weight b there are N_(G-k)(b) possible nonzero
patterns. The selected fields supply at most N_k(d_s-b)
choices. Sum that many largest upper terms within the
exact mask capacities, multiply by N_(G-k)(b), and sum over
b for which those counts can be positive. Multiply by
binom(G,k) for all chosen subsets. Sum over k=1..G and
d_s=G..2t to obtain the integer correction R_t.

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
resulting first-pair errors. There are four cell factors, two each of widths 24 and 1. Each possible
systematic error determines its direction uniquely. The remaining
a-1 systematic and a parity fields have u=(2a-1)*50 bits.
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
this field. The coefficient comes from requiring all k mask
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
product of their full weight enumerators, and denote its
weight-w coefficient by H_w. The product histogram gives
the number of tuples with each polynomial H. Under fixed
local total noise weight w and that syndrome tuple, a
true first-pair weight d occurs in exactly H_d*binom(u,w-d)
true noise patterns. Let V_H(w) be the sum of these counts
over d. The orthogonal coordinates outside this run do
not enter the local comparison.

On a distance tie, swapping the true word and its run
competitor at the same received word preserves the run,
the local total noise weight, the unchanged noise outside
the run, and the full first-pair syndrome tuple. At the
first pair the cube predecessor is unchanged, and local
systematic and parity differences have the form e,L(e).
This swap is a fixed-point-free involution on ordered
tied competitors and reverses the lexicographic order.
Within each class, exactly half of the tied competitors
are oriented toward the earlier word. Strictly nearer
competitors contribute their full counts. A union bound
over the resulting oriented pairs upper bounds the number
of removed true pairs. True noise patterns in the class
have equal probability conditional on their local weight,
and actual messages are uniform.

Given a specific true first-pair pattern of weight d0,
there are H_d candidate first-pair patterns of weight d
except that the one unchanged direction is omitted at
d=d0. Every other first direction satisfies the uniform
rank bound above. Therefore the doubled denominator's
numerator for candidate counts averaged over messages
and summed over the true noise patterns in this class
is at most

    V_H(w) * sum_(d<=w) H_d * D_(a,w-d)
        - sum_(d<=w) H_d * binom(u,w-d) * D_(a,w-d).

The subtraction omits the specific first direction equal
to the true word before applying the uniform bound, one
such possibility per actual pattern. Clip this integer
at 2^(p+1)*V_H(w), since the number of removed true
pairs per message is at most the number of true patterns.
Add the clipped bounds over H with the profile-product
multiplicities. Divide by 2^(p+1) to get the upper
count of removed true patterns of local weight w, averaged
over messages. Multiply this count by the number of
outside noise patterns at weight t-w and sum over w.
This is the average upper count for a fixed start and run
length in total shell t. All five start positions have the
same count by cyclic symmetry. Sum over the four possible
lengths and the five starts before rounding to an integer.
All operations are exact integer or rational arithmetic.

## Occupancy after removals by exact noise weight

At each total noise weight v, upper-bound the mean number
of removed patterns per true message. Remove a pair of true
and received words if the anchor or child admission test fails,
or if some single-run competitor is closer, or ties and is lexicographically earlier. Two retained
words that differ on a proper subset cannot coexist, since the
farther is removed, or the later one on a tie.

For the anchor, sum

    binom(ANCHOR_BITS,d)*binom(500-ANCHOR_BITS,v-d)

over d above the radius and at most v. The child count was
computed for the disjoint range within the radius. For the
proper-run term, use the convolution over unaffected noise
coordinates described in the previous section. Before
multiplying by the common channel probability
9^(500-v)/10^500 its value is the upper mean count of
removed patterns at weight v. Sum over the run lengths and
their starting positions, then round upward to an integer.
These terms all bound mean counts under the uniform message.
Add the three upper counts and clip at binom(500,v), the
whole shell. Denote this integer by U_v. This upper bound
allows repeated counting when several reasons remove a pair.

At cutoff t<=65 let V_t be the volume of the 500-bit ball,
and A_t=V_t-sum_(v<=t) U_v. At least Q*A_t ordered pairs of
true and received words remain, for Q=2^250. For each
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
9^(500-t)/(Q*10^500). Summation by parts with the
nonincreasing weights 9^(500-t), and ignoring successful
decodes beyond the cutoff, gives success at least

    sum_(t=0..65) (G_t-G_(t-1))*9^(500-t)/10^500.

The four components use disjoint noise and uniform
independent messages. The product decoder succeeds when
all components succeed. Thus its error is at most one
minus the fourth power of the component success bound.

Part 3 checks the claimed error bound from this exact lower
success bound and the validated witness counts.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
For each component the ten received fields are copied to padded little-endian
one-word fields in arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
Words SAVED_B through SAVED_B+14 retain the first three component results.
The remaining working arrays are consecutive and follow these saved words.
The fixed data include reduction and square tables, the repetition mapping and candidate tables, and four additional
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
- BUDGET_BASE, one word per anchor depth 0 through 11 holding C(d).
- PROFILE_BASE[m], one word per orbit (width 24) or syndrome
  (width 1) giving its complete-profile type number.
- CAP_OFFSETS[i-1], one block per emitted cell i, with entries
  E_(i,p)(r) at offset p*(1051)+r for every type p and state r.

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
can be discarded. All 16 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use two accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=1 word and h=50 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^50. This fits in one word since the
largest possible degree is 12. To multiply a field
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
c*X^50 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=7 byte slots, emit a
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
Move R registers to zero. For each of n0=13 nibble
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
For component k=0, 1, 2 or 3 emit a literal copy from
systematic starting bit 250*k and parity starting bit 1000+250*k.
Each group contains five fields of 50 bits. Use the extraction
routine for the ten padded words of Y and P, then store each.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 30. Store zero in the five
padded incumbent words. Copy Y0 to X0 in two primitives,
set Depth to zero and store NEXT[0]=0. This setup costs at most
108 primitives per component. Each visited cell and anchor frame
writes every slot it later reads.

For packing a field from known short pieces, clear one
register by a free move and insert each piece with a left shift
and OR at its literal bit offset. Store the accumulator at
the field's word address.

Emit the complete component traversal below four times, with
the corresponding literal input offsets and fresh setup.
After each of the first three returns, copy the five padded
incumbent words B to SAVED_B+5*k using ten primitives.
After the last return the sixteen output words gather the first
750 bits from SAVED_B and the last 250 bits from the current B.
Each output word receives pieces from at most three successive
padded source words, since the source fields have width 50.
Extract each literal piece and insert into a zeroed output
accumulator; this costs at most 3*(6+2)+1 including the output
store. The last output word has its high 24 bits zero. Combined
non-traversal work including setup, transfers and stores
is at most 4*108+3*10+16*25 <1000.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=50, and Depth is a persistent register.

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
        branch if done to this component’s return continuation
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
N=sum_(i<=R) binom(50,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base (QR_PARITY_BASE at width 24),
load the parity image, and pack
these images into a one-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining by adding Depth to the literal BUDGET_BASE and loading
that word. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 250 further primitives:
the four cells use at most 4*(6+2+4), and finishing AP
uses 3, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 16 cells in field and within-field order.
Its frame word slots at offsets 0..5 are Base, U, Rem,
Dist, Rank and Perm. The literal SELECTED slot holds
this cell's selected systematic value. Parent continuation
means the previous block's CELL_NEXT, or ANCHOR_EXTEND
for the first block.

For width m=1 use the following emitted program.

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
        j := mul(p,1051)
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
        needed := add(cc,16-DepthCell)
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
        run FIELD after each of the three nonfinal fields
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
    j := mul(p,1051)
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

At FIELD and FINAL, pack the four selected cell slices from
the current completed field into its one padded Xj word,
with at most 4*(1+4)+1=21 primitives. Compute Cj=Xj^3.
At FIELD this cube supplies the next field; its total routine
charge is at most CubeOps+180. At FINAL add the distance of
the anchor parity to Distance by loading AP and C4 once, XORing, popcounting and adding.
The residual AP+C4 is P0+L(X0)+C4. Distance already includes
the anchor systematic distance and all four visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all five systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 200 primitives.
## Exact worst-case operation bound

Number the 16 literal cells by i=1..16. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 67 at width 24 and 47 otherwise, and g_i be 80 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
T_i is CubeOps+180 at i=4,8,12 and CubeOps+200 at i=16.
It is zero at all other indices.

For all possible states r=0..1050 define in decreasing cell order

    F_17(r)=0
    U_(i,p)(s)=sum_(c<=s-(16-i))
                 A_(p,c) * (g_i + T_i + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,1050]: U_(i,p)(s)<=F_i(r)-h_i}.

The nonnegative recurrence is nondecreasing in r by induction
and U is nondecreasing in s. The maximum defining E exists
because r itself is allowed. Its search returns the rightmost
valid index, including equality. These are the fixed CAP table entries.
For any actual state r and profile p, the accepted attempts
are exactly the rows up through their recorded cost threshold
E_(i,p)(r)-(16-i). Their continuation states E_(i,p)(r)-c
belong to [0,1050] because costs are positive and successful
attempts satisfy the reserve. At a successful attempt the
remaining program cost is bounded by g_i+T_i+F_(i+1) at that state.
Consequently the actual entire cell traversal costs at most
h_i+U_(i,p)(E_(i,p)(r))<=F_i(r). This proves the bound by
downward induction without distributional hypotheses.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 250+CubeOps. Hence a worst-case
all-input integer bound for all four traversals and their layout is

    W = 1000 + 4 * sum_(d=0..11) binom(50,d) *
                 (44 + 250 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.

## Scores

The exact counted primitive bound is W=434230519218960832, with
W^125 <= 2^7324. The verifier checks this
inequality and the exclusive RAM address limit.

The exact success calculation above gives error at most
B=728158332970513879/1000000000000000000. The declared experiment checks that rational inequality
using integer or fraction arithmetic throughout.
