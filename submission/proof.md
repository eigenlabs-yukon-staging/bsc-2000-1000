# Two cyclic cube fields with weight-conditioned admission

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

The decoder enumerates all anchors x0 of distance d<=98
from received Y0. For each use C(d)=4635-40d as the integer
child cost allowance. Starting from the first child cell,
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
use a row of the first 2^23 noise patterns, through an explicit
syndrome orbit map and coordinate permutation defined below.
This preserves the distances and the syndrome. The selected
cell is u+e, where e is the systematic part of its noise pattern
in original coordinates. All rows are in nondecreasing distance
and thus nondecreasing cost.

Depth-first search the 23 cells, retaining a prefix of length i
only at cost at most C(d)-(23-i). Since future costs are at
least one, every full tuple of cost at most C(d) using the fixed
rows passes all tests. Score each complete tuple by full
received distance. Use the all-zero fallback and only replace
the incumbent on a strict improvement, with the traversal order
below. The error proof and the work count use exactly this rule.

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

## Exact child-cost distributions

On the true path the predecessor cube equals the true cube and
subtracting it leaves the original independent observations.
For the channel with crossover probability 1/q, where q is one
of 10, 12 and 14, a particular cell noise of weight w in width m
has probability (q-1)^(2m-w)/q^(2m). The decoder always uses
the nominal cost defined from the weight base 9. The auxiliary
channels only bound admission probabilities.

For each cell width and syndrome, take the same actual retained
pattern count as the fixed decoder row. At width 24 this count
can be truncated inside the distance tie at the last permitted
rank. At widths 4 and 12 it is the full row enumerator. With that
count, add its product with (q-1)^(2m-w) at the assigned nominal
cost. Summing over the local syndromes and true weights gives
the numerator of the true-cost law, with missing mass for the
unlisted rows. Multiply twenty width-24, one width-12 and two
width-4 factors and truncate above 4635. The denominator is
q^1000. Coefficients are nonnegative and bounded by the total
numerator at each step, at most q^1000, so the bit length of
q^1000 (rounded up to a byte boundary) suffices per digit in the
packed integer calculation.

For anchor systematic distance d<=98 define E_d to be the
event that the correct child path is absent from its rows or
exceeds total cost C(d). By the prefix reserve test these are
exactly the reasons it fails. Its probability under the 1/q
child noise is one minus the appropriate coefficient prefix
divided by q^1000. The arrays in tails.json contain these
probabilities rounded up to denominator 10^60, for all three
channels and all depths 0..98. Declared parts 1 and 2 check
every entry. The event depends only on child channel noise, not
on either anchor noise or the actual message, because on the
correct path the adjustment and the linear local parity give
the true noise syndrome and the true error pattern.

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
Call the resulting integer D_t. The expected
number of full-change competitors within total
t for each fixed true noise is at most D_t/Q.
The checker verifies D_t is nondecreasing through t=230.

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

## Occupancy of the actual candidate list

At distance cutoff t<=230 consider all pairs of true word and
received word with true distance at most t. Remove any pair
for which the true word fails the specified anchor or child
admission tests, and any pair for which a single-run competitor
has distance no greater. The fallback is allowed even when it
would have been removed. Let L_t count removed pairs, Q=2^1000,
and V_t be the 2000-bit Hamming ball volume through t. Then
Q*V_t-L_t pairs remain. Two distinct retained candidates at the
same received word must differ in both fields, since otherwise
the farther would be removed, or both would be removed on a tie.
Their number of ordered distinct pairs at common received words
is at most V_t*D_t by summing the per-noise full-change bound.

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
success is 9^(2000-t)/(Q*10^2000). Thus the success through the
cutoff is obtained by summation by parts with w_t=9^(2000-t)
and w_231=0, by multiplying F_t/Q by (w_t-w_(t+1))/10^2000
and summing. Ignoring success outside the cutoff only decreases
the lower bound. A removed pair with true noise weight v can
be charged only at t with i(t)>=v, and alpha_(i(t))<=alpha_v
because D_t is nondecreasing. Its total charge is at most
alpha_v*w_v/(Q*10^2000) by telescoping.

The error bound therefore consists of the unremoved baseline

    1 - sum_(t<=230) (U_t-U_(t-1))*w_t/10^2000,

plus two proper-run terms (one per start), plus the admission
failure charge weighted by alpha at the true total noise weight.
Here U_(-1)=0. The proper-run terms use the J_(1,v) expression
above, multiplied by alpha_v and summed through 230. An upper
bound on the weighted admission charge may count both removal
reasons even when they occur on the same pair.

## Weight-conditioned admission charge

Write d for the true anchor systematic weight on 500 bits.
Write w for its independent parity noise on the other 500 anchor
bits and Z for the child total noise on 1000 bits. The total weight
is d+w+Z, with the three noise blocks mutually independent.
If d>98 the anchor is absent whatever the other noises; its
exact weighted charge is

    sum_(d>98) Pr[Bin(500,1/10)=d]
        * sum_(z<=230-d) Pr[Bin(1500,1/10)=z] * alpha_(d+z).

Use only d<=230 because the weighted charge vanishes otherwise.

If d<=98, failed admission is exactly E_d under the child's
nominal noise. For q=10,12,14 the likelihood ratio of a child
noise pattern with weight z under 1/10 to its probability under
1/q is

    R_q(z)=q^1000*9^(1000-z)/(10^1000*(q-1)^(1000-z)).

It is nondecreasing in z, since its step multiplier is (q-1)/9.
Hence for arbitrary E_d,

    Pr_(1/10)[E_d and Z<=z] <= min_q R_q(z)*Pr_(1/q)[E_d].

This inequality covers both failure to find a required QR row entry
and failure of the cost test. It does not assume a monotone relation
between noise weight and normalized cost.

Put alpha_231=0 and define

    G(s)=sum_(w<=230-s) Pr[Bin(500,1/10)=w]
                            *(alpha_(s+w)-alpha_(s+w+1)).

For fixed d,w, summation by parts over z turns the truncated
weighted child charge into nonnegative multiples of the cumulative
probabilities Pr[E_d and Z<=z]. Independence of w permits summing
first over it, producing G(d+z). The child contribution is thus at
most the sum for d<=98 and z<=230-d of

    Pr[Bin(500,1/10)=d] * G(d+z) * min_q R_q(z)*Pr_(1/q)[E_d].

The verifier computes the exact anchor term. For the child term
it rounds G up to denominator 10^60, multiplies the verified tail
ceilings by the exact ratios, and rounds these products up to
denominator 10^60 before minimizing. Each replacement is upward,
and both weights are nonnegative. It then sums the resulting
integer numerators against the exact anchor binomial masses.
Part 3 checks this weighted error bound together with the full
change and proper-run contributions.

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
and Hamming mapping and candidate tables, and four additional
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

Construction.py and qr.py specify every fixed datum. The
exclusive end of RAM is 3703713469, within the model's
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
word, XOR and store to AP. Set Remaining=C(Depth) using one multiplication by 40 and one
subtraction from 4635. Move Distance=Depth and jump to the first
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
        store(Rem,Remaining)
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
Initialization uses at most 29. The unique stopping attempt
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
    pos := mul(perm,6*256)
    pb := add(QR_PERM_BASE,pos)
    store(Base,base)
    store(Perm,pb)
    store(U,u)
    store(Rem,Remaining)
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
The QR initialization uses at most 36 primitives, including
its three initial extracts, and its single failed attempt at
most 11. Charge at most 60 for these together and at most
80 for each successful attempt including its jump. There is
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
If replacing, move BestD to TotalD and copy all sixteen
padded systematic words of X to B. This is at most 320
primitives plus the cube. It retains the earlier candidate
on equality and writes the 1000 systematic bits as prescribed.

## Exact worst-case operation bound

Number the 23 literal cells by i=1..23. Let A_(s,c) count
the candidates of cost c in each possible row of the cell's
width, from the checked full syndrome histogram and the
same row cap as the fixed decoder table. Since the width-24
permutation preserves distance, its counts are exactly the
truncated histogram counts. At any actual parent let r be
the remaining budget. Its accepted costs satisfy c<=r-(23-i).
The checker computes the bound recurrence

    F_24(r) = 0
    F_i(r) = h_i + max_s sum_(c<=r-(23-i))
                          A_(s,c) * (g_i + T_i + F_(i+1)(r-c)).

Here (h_i,g_i)=(60,80) for width 24 and (40,28) otherwise.
T_i is zero before the final cell, CubeOps+320 at the
final cell. By downward induction F_i(r) bounds all primitives
from CELL_INIT through final return to its parent. Its
fixed charge counts setup and the unique stopping attempt.
Each accepted attempt adds the cost of the dependent child
with exactly its new budget. Every actual syndrome is included
in the enumerated maximum. Allowing further independent
maximization at downstream nodes is valid whatever candidates
were selected earlier. The explicit field placement and jump
charges already appear in the local cost constants.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 560+CubeOps. Hence a worst-case
all-input integer bound is

    W = 1000 + sum_(d=0..98) binom(500,d) *
                 (44 + 560 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.


The declared exact output gives the following bounds, with the
probability components rounded upward at denominator 10^18.

    Unremoved baseline error <= 55175259434787520/10^18
    Single-start proper-run charge <= 1/10^18
    Anchor admission charge <= 329897/10^18
    Child admission charge <= 498567716/10^18
    B = 13793814983421283/250000000000000000
    W = 1754609142709721328150812810937476980479266661579811597200167634668164187114807275581680758013414239156603791189562581960857276
    W^(8) <= 2^(3355).

The exact sum is rounded upward after combining the unrounded
components. Independent rounded component displays may sum to
slightly more than the final combined rounded bound.
