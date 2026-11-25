# Six cyclic cube fields with two shortened coordinates

## Encoder and restricted decoder

Work in F2[X]/(P) with P=X^167+X^6+1 and the polynomial basis,
low coefficient first. Degree 167 is prime, so Rabin uses
X^(2^167)=X modulo P and gcd(P,X^2+X)=1. Let U consist of the field
values whose two highest coefficients are zero.

Split the 1000 message bits into a 165-bit field value b0 in U
and five 167-bit field values b1 through b5. The six parity fields use

    p_i = L(b_i) + (b_(i-1))^3.

Parity field p1 has its two highest bits omitted;
all others transmit all 167 bits. Thus there are 1000 systematic
and 1000 parity channel bits. The encoder is injective by its
systematic prefix. Construction.py fixes all layouts.

Each L consists of six width-24, one width-12, two
width-4, and three width-1 identity cells in increasing bit order.
For a width-24 input u use low 23 parity bits equal to the
remainder of u*X^23 modulo the degree-23 binary G=9205487.
Set its high parity bit to the parity of the weight of u and
its remainder. The width-12 Golay cell uses the same rule
with degree 11 and polynomial 1+X^2+X^4+X^5+X^6+X^10+X^11.
The width-4 cell has columns 14,13,11,7 and the three
width-1 parity cells are identities. All maps are linear.
In particular the high two coordinates of L(b0) are zero.

The decoder visits each anchor x0 in U within 34 systematic
errors of Y0. At depth d (the anchor systematic distance), use
C(d)=2869-48d. For a fixed anchor, visit child fields 1 through 5
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

There are 60 total child cells. Continue with a prefix of length
i only if its used cost is at most C(d)-(60-i). Every complete
tuple of cost at most C(d) passes all prefix tests, since every
later cell costs at least one. Compare all returned messages by
their full distance. Initialize the incumbent to zero, replace
only on a strict improvement, and use the traversal order below.
This defines the selected message on every received word.

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
of 10 and 12, a particular cell noise of weight w in width m
has probability (q-1)^(2m-w)/q^(2m). The decoder always uses
the nominal cost defined from the weight base 9. The auxiliary
channels only bound admission probabilities.

For each cell width and syndrome, take the same actual retained
pattern count as the fixed decoder row. At width 24 this count
can be truncated inside the distance tie at the last permitted
rank. At the other paired widths it is the full row enumerator. With that
count, add its product with (q-1)^(2m-w) at the assigned nominal
cost. Summing over the local syndromes and true weights gives
the numerator of the true-cost law, with missing mass for the
unlisted rows. For each of the two unpaired bits include a one-bit factor with
numerator q-1 at cost c(10,9) and numerator 1 at c(10,1).
Multiply thirty width-24, five width-12, ten width-4, and
thirteen paired width-1 factors and the two unpaired factors. Truncate
above 2869. Their shared denominator is q^1668. Every product
has coefficient sum at most q^1668, so its bit length (rounded
up to a byte boundary) suffices per packed coefficient.

For anchor systematic distance d<=34 define E_d to be the
event that the correct child path is absent from its rows or
exceeds total cost C(d). By the prefix reserve test these are
exactly the reasons it fails. Its probability under the 1/q
child noise is one minus the appropriate coefficient prefix
divided by q^1668. The arrays in tails.json contain these
probabilities rounded up to denominator 10^60, for both channels and all depths 0..34. Declared parts 1 and 2 check
every entry. The event depends only on child channel noise, not
on either anchor noise or the actual message, because on the
correct path the adjustment and the linear local parity give
the true noise syndrome and the true error pattern.

## Competitors changing all six fields

For a full degree-167 field, write q=2^167. For e nonzero,
the cube difference is

    (x+e)^3+x^3 = e^3 + e*x^2 + e^2*x.

The map z -> z^2+z has kernel {0,1} and image the kernel
of the absolute trace. The trace is an F2-linear nonzero
functional by the polynomial root bound. Thus for x uniform
over the full field, the parity difference is uniform on an
affine hyperplane. Its normal has entries Tr(e^(-3)*X^j).
The pairing (h,z) -> Tr(hz) is nondegenerate because a nonzero
h permutes the field and trace is nonzero. Thus a normal
determines the field multiplier e^(-3). Since 167 is odd, cubing permutes the field and there is
at most one direction per normal.

For parity p1 the predecessor x lies in U, and only the low
165 coordinates of its difference are transmitted. Any mask
on the transmitted coordinates lifts to a unique multiplier h
for the full trace pairing, by extension as a functional
on the 167 basis coordinates. For a nonzero mask h is nonzero.
It is constant on the derivative image of U exactly when

    L_e^*(h) = (h*e)^(2^166) + h*e^2

belongs to U-perp, the orthogonal complement under trace.
Indeed the displayed adjoint follows because
Tr(h*e*x^2)=Tr((h*e)^(2^166)*x). The space U-perp has
four elements. Fix any z from it and square the equation
L_e^*(h)=z to obtain

    h*e + h^2*e^4 = z^2.

This is a nonzero degree-four polynomial in e, with at most
four roots. For each nonzero mask there are therefore at most
16 possible nonzero directions for which it occurs. No
assumption about the actual projection rank is needed for this
per-mask bound.

When all fields change, fix the differences e_i. Parity
differences in the six transmitted fields are independent because
they use distinct true predecessor fields, which are independent
under the uniform systematic prior. Within each field the
distribution is uniform on a translate of a binary linear image.
If that field has parity length n, its probability at a vector
is 2^-n times a sum of characters over all annihilating masks,
with translation signs. For example on an affine support t+V
the sum of the characters indexed by V-perp is |V-perp| at x
in the support and zero outside; dividing by 2^n gives 1/|V|.
The linear term L(e_i) is just a translation, including after
puncturing. Multiplying all six expansions gives a factor 1/Q
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

Suppose the term selects k of the five ordinary parity masks.
Optionally it also selects a nonzero mask on p1, indicated by
h=1 (otherwise h=0). There are binom(5,k) choices of ordinary
positions. A combined weight-w mask appears in at most the
following counts across the corresponding predecessor directions:

    h=0: [Z^w] ((1+Z)^167-1)^k
    h=1: 16*[Z^w] ((1+Z)^167-1)^k * ((1+Z)^165-1).

These follow from ignoring zero masks at selected positions
and the one or sixteen direction preimages. The calculation
evaluates those coefficients by inclusion and exclusion.
There are n=167k+165h systematic bits at the selected
predecessor positions, with the other 1000-n systematic bits
unrestricted by these masks. For fixed true noise and
combined candidate systematic error weight a on selected
predecessors, there are at most binom(n,a) choices of their
directions, since each possible systematic error determines them.
For a fixed nonzero direction e in U the linear derivative on U
has kernel {0,e} and image dimension 164. Projection can
remove at most two dimensions, so its rank is at least 162
on the 165 parity coordinates. There are consequently at most
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
1000 parity bits. Call the resulting integer D_t. The expected
number of full-change competitors at distance at most t, for every
fixed actual noise, is bounded by D_t/Q. The checker verifies
the bound's nondecreasing values through t=230.

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
therefore have twelve small factors in the first-pair polynomial H.

Let p be the total length of parity fields i+1 through i+a,
let u be p plus the systematic length of the other changed fields,
and let n be u plus both lengths in the first pair. The parity
differences among these p bits are independent for fixed directions.
If no changed predecessor is b0, their point mass is at most
2^a/2^p. If b0 changes, its direction e is a nonzero element
of U. The map e*x^2+e^2*x has kernel {0,e} on U, so
its image on the 165-dimensional domain has dimension 164.
Puncturing can lose at most two independent coordinates and
leaves image dimension at least 162. Its point mass on the
165 transmitted coordinates is at most 8/2^165. Thus in either
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

For each geometry and run length, the fixed GRID in verify.py
lists rational pairs (r,j). Set z=(r/100)^8 and rho=j/8. The
coefficient prefix is at most z^-w*(1+z)^u*H(z). The chance
of any competitor is bounded by min(1,x) for the expected
count x, and min(1,x)<=x^rho for 0<rho<=1. Thus use

    2^(-(p-a-2*k)*j/8) * (r/100)^(-j*w) * (1+z)^(u*j/8)
              * product_cells H_s(z)^(j/8),

where k indicates that b0 is among the run's changed fields.
For a nonnegative rational A/B, power_up computes the upward
integer approximation ceil(2^64*(A/B)^(j/8)). It obtains the
floor by three integer square roots of floor(2^512*A^j/B^j).
Compare the eighth power times B^j against 2^512*A^j and
increment unless equal. Use this bound separately for each
local H_s(z) factor. Also compute the factors 2^(-(p-a-2*k)*j/8)
and (1+z)^(u*j/8) by splitting their exponents into the
integer part and a remainder over 8, using the exact
same upward eighth-root operation for the fractional part.
The factor (r/100)^(-j*w) is exact rational.

For each first-pair cell multiply the rounded root numerator
by h_s(d) and sum over the syndromes at each true weight d.
This covers all actual local noise patterns, including the
fixed-parity and unpaired systematic cases above.

Form the product of the twelve polynomials so obtained and
(1+Z)^u. At each local weight w multiply by the rational
common factors above, dividing by the twelve powers of 2^64
for the local roots and the two additional powers of 2^64
for the common rounded factors. This bounds the sum of the individual probability
bounds over true local noise patterns of that weight. Minimize
the result at each w over the stated grid and the trivial
binom(n,w) bound. Each candidate is uniform for every fixed
outside noise of that weight class, so the minimum before outside
expansion is valid. The local sum is then multiplied by
binom(2000-n,v-w) and channel factor 9^(2000-v)/10^2000.
Sum over w to get a bound J_(i,a,v) on each joint probability
for start i, length a and total true weight v. The checker
groups exactly identical geometries by their multiplicities.

The root numerator is at most 2^88 since each cell has no
more than 2^24 terms at z<=1 and j<=8. Each product of paired
components sums to at most 2^(334+12*88); the free factor
adds at most 2^1503. Thus 4000-bit digits suffice throughout
the packed polynomial products. All further arithmetic is rational.

## Occupancy of the actual candidate list

At distance cutoff t<=230 consider all pairs of true word and
received word with true distance at most t. Remove any pair
for which the true word fails the specified anchor or child
admission tests, and any pair for which a single-run competitor
has distance no greater. The fallback is allowed even when it
would have been removed. Let L_t count removed pairs, Q=2^1000,
and V_t be the 2000-bit Hamming ball volume through t. Then
Q*V_t-L_t pairs remain. Two distinct retained candidates at the
same received word must differ in all six fields, since otherwise
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

plus the proper-run terms for all five lengths and six starts, plus the admission
failure charge weighted by alpha at the true total noise weight.
Here U_(-1)=0. The proper-run terms use the J_(i,a,v) expression above,
weighted by alpha_v and summed through 230 for each start and length. An upper
bound on the weighted admission charge may count both removal
reasons even when they occur on the same pair.

## Weight-conditioned admission charge

Write d for the true anchor systematic weight on 165 bits.
Write w for its independent parity noise on the 167 anchor
bits and Z for the child total noise on 1668 bits. The total weight
is d+w+Z, with the three noise blocks mutually independent.
If d>34 the anchor is absent whatever the other noises; its
exact weighted charge is

    sum_(d>34) Pr[Bin(165,1/10)=d]
        * sum_(z<=230-d) Pr[Bin(1835,1/10)=z] * alpha_(d+z).

Use only d<=230 because the weighted charge vanishes otherwise.

If d<=34, failed admission is exactly E_d under the child's
nominal noise. For q=10,12 the likelihood ratio of a child
noise pattern with weight z under 1/10 to its probability under
1/q is

    R_q(z)=q^1668*9^(1668-z)/(10^1668*(q-1)^(1668-z)).

It is nondecreasing in z, since its step multiplier is (q-1)/9.
Hence for arbitrary E_d,

    Pr_(1/10)[E_d and Z<=z] <= min_q R_q(z)*Pr_(1/q)[E_d].

This inequality covers both failure to find a required QR row entry
and failure of the cost test. It does not assume a monotone relation
between noise weight and normalized cost.

Put alpha_231=0 and define

    G(s)=sum_(w<=230-s) Pr[Bin(167,1/10)=w]
                            *(alpha_(s+w)-alpha_(s+w+1)).

For fixed d,w, summation by parts over z turns the truncated
weighted child charge into nonnegative multiples of the cumulative
probabilities Pr[E_d and Z<=z]. Independence of w permits summing
first over it, producing G(d+z). The child contribution is thus at
most the sum for d<=34 and z<=230-d of

    Pr[Bin(165,1/10)=d] * G(d+z) * min_q R_q(z)*Pr_(1/q)[E_d].

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
The twelve received fields are copied to padded little-endian
three-word arrays Y and P. Candidate fields X, their cubes C, and
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

Construction.py and qr.py specify every fixed datum. The exclusive end of RAM is the value MEMORY_END in
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
can be discarded. All 60 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use three accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=3 words and h=39 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^167. This fits in one word since the
largest possible degree is 9. To multiply a field
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
c*X^167 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=21 byte slots, emit a
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
Move R registers to zero. For each of n0=42 nibble
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
Use this for the 36 padded words of Y and P, using the
transmitted widths in construction.py, then store each.
The shortened Y0 and P1 have two extra zero padding bits.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 108. Store zero in the 18
padded incumbent words. Copy Y0 to X0 in six primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear three
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the three
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. For output use field widths 165,167,167,167,167,167. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=165, and Depth is a persistent register.

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
N=sum_(i<=R) binom(165,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base (QR_PARITY_BASE at width 24),
load the parity image, and pack
these images into a three-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining=C(Depth) using one multiplication by 48 and one
subtraction from 2869. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 250 further primitives:
the 12 cells use at most 12*(6+2+4), and finishing AP
uses 9, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the 60 cells in field and within-field order.
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
        needed := add(cc,60-DepthCell)
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
        run FIELD after completing each of the four nonfinal fields
        run FINAL if this is the literal final cell
        jump to the next CELL_INIT or, at the final cell,
            back to its CELL_NEXT

Each nested load or arithmetic expression is a separate primitive.
Initialization uses at most 29. The unique stopping attempt
at a visited frame costs at most 11. Each successful attempt
costs at most 28 including its jump, plus FIELD or FINAL where applicable.

For the two systematic-only cells of field 1, use width one.
Instead of the ordinary initialization, extract the systematic
bit and perform the same five stores, using the literal
SYS_ONLY_BASE as Base. The row has exactly two candidates,
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
80 for each successful attempt including its jump. Each field boundary falls at an width-1 cell. All stages use the
six frame words and fewer than 20 scratch registers alongside
the persistent registers specified above.

When an accepted child returns, the parent still has its
Base, U, saved budget and distance and next rank in RAM.
The frame stops when it runs out of entries or finds the
first over-budget record. The costs are sorted, so further
entries would also fail. Each full tuple in the specified
list is therefore visited exactly once.

At FIELD and FINAL, pack the twelve selected cell slices from
the current completed field into its three padded Xj words,
with at most 12*(1+4)+3=63 primitives. Compute Cj=Xj^3.
At FIELD this cube supplies the next field; its total routine
charge is at most CubeOps+180. At FINAL add the distance of
the anchor parity to Distance by loading AP and C5 round by
round over three words, XORing, popcounting and adding.
The residual AP+C5 is P0+L(X0)+C5. Distance already includes
the anchor systematic distance and all five visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all eighteen systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 200 primitives.
## Exact worst-case operation bound

Number the 60 literal cells by i=1..60. Let A_(s,c) count
the candidates of cost c in each possible row of the cell's
width, from the checked full syndrome histogram and the
same row cap as the fixed decoder table. Since the width-24
permutation preserves distance, its counts are exactly the
truncated histogram counts. At any actual parent let r be
the remaining budget. Its accepted costs satisfy c<=r-(60-i).
The checker computes the bound recurrence

    F_61(r) = 0
    F_i(r) = h_i + max_s sum_(c<=r-(60-i))
                          A_(s,c) * (g_i + T_i + F_(i+1)(r-c)).

Here (h_i,g_i)=(60,80) for width 24 and (40,28) otherwise.
T_i is CubeOps+180 at i=12,24,36,48 and CubeOps+200 at i=60.
It is zero at all other indices. By downward induction F_i(r) bounds all primitives
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
anchor cube cost at most 250+CubeOps. Hence a worst-case
all-input integer bound is

    W = 1000 + sum_(d=0..34) binom(165,d) *
                 (44 + 250 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.






The declared computation gives the following ceilings with the
common denominator 10^18 for the probability components.

    Unremoved baseline error <= 55793697426781823/10^18
    Length-wise proper-run charges over all starts <= (2670372810356, 331045183377248, 930975720879573, 556912090598107, 27568623152171)/10^18
    Anchor admission charge <= 8964060087830/10^18
    Child admission charge <= 2976690510025831/10^18
    B = 12125704797542587/200000000000000000
    W = 363795158065148346627881170867534718892059028777348856557051491376
    W^(1000) <= 2^(217789).
    Exclusive RAM end = 3703600094.

For B, round the exact sum upward after combining the unrounded
terms, whose length charges already include every start. Separately rounded
component displays may sum slightly above the combined ceiling.
