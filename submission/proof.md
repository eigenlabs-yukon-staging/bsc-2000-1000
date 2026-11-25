# A mixed product of two-field and three-field cyclic components

## Encoder and decoder

Use the polynomial basis of F2[X]/(P), where
P=X^50+X^9+X^2+X+1. The declared field check applies
the Rabin criterion with prime divisors 2 and 5 of 50.

Split the 1000 payload bits into six 150-bit components
followed by a 100-bit component. In a component write g
for its number of 50-bit fields, K=50*g for its payload,
and n=2*K for its code length. The first six components
have g=3 and the last has g=2. Within a component its
systematic fields b_i have increasing little-endian order.
The parity fields are p_i=L(b_i)+b_(i-1)^3, with indices
cyclic modulo g. Transmit all systematic fields in component
order, then all parity fields in the same order. This is
a systematic injection on 1000 bits into 2000 bits.

For each component put J=4*(g-1) for its number of child
cells. The two entries THREE and TWO in construction.py
give the literal anchor budget list C(d), anchor radius R
equal to the length of that list minus one, and distance
cutoff T. The cutoffs are 48 for g=3 and 38 for g=2.
The cutoff is used in the error calculation. The decoder
scores all complete tuples admitted by its actual search.

In each field L is block diagonal with two 24-bit cells and two
width-1 repetition cells in increasing bit order. For a width-24
input u, the low 23 parity bits are the binary polynomial
remainder of u*X^23 modulo the degree-23 G=9205487 (in
low-coefficient-first order). The high parity bit is the parity
of the sum of Hamming weights of u and the low parity bits.
The width-1 parity map is the identity. These maps are linear. The file construction.py fixes the integer layout.

The decoder enumerates all anchors x0 of distance d<=R
from received Y0. For each use entry C(d) at index d of the component budget list as the initial
remaining-cost state. Visit the child fields j=1 through g-1 in order. In each field,
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

Depth-first search the J child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(J-i), where r is
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

## Exact admitted child noise by weight

On the true path the predecessor cube equals the true cube.
Subtracting it leaves the original independent observations.
For any local syndrome of a given profile type p, the number
of retained row entries of true weight w is determined by its
profile and row cap. The orbit maps preserve weights. If the
cap cuts a distance tie, its retained count at that weight is
still the prescribed count. The type multiplicity is the
number of distinct syndromes with that profile.

For each cell position i, type p, nominal cost c and local
weight w, let H_(i,p,c,w) be the retained count at that cost
and weight, multiplied by the number of syndromes of the
type. In a fixed cell the mapping between its raw noise
pattern and its syndrome together with its systematic noise
is bijective. Thus H counts all raw local noise patterns that
would select one of the designated entries on the true path.
Compute it from the exact verified profile and the integer
nominal cost. For a given state r and type, the effective
cap e=E_(i,p)(r) is the same as in the RAM schedule and the
entry is accepted precisely if c<=e-(J-i). Children on
the true path depend only on their own noise, so a disjoint
product of local choices counts all admitted sequences.

Give every accepted suffix a monomial whose degree is its
total noise weight. With the empty suffix P_(J+1,r)(X)=1,
the exact recurrence is

    P_(i,r)(X) =
        sum_(p,c,w with c<=E_(i,p)(r)-(J-i))
        H_(i,p,c,w) X^w P_(i+1,E_(i,p)(r)-c)(X).

For anchor depth d<=R the coefficient of X^z in
P_(1,C(d)) is the number A_(d,z) of admitted child noise
patterns of weight z on n-100 bits. Exponents and all weights
are nonnegative. Therefore coefficients above the component
cutoff T can be discarded after each step. The verifier
packs one degree per digit of width n-99 bits, since all
coefficients even before truncation are bounded by the
number 2^(n-100) of raw child noise patterns. No product
or sum can carry into another digit. Shifting by w digits
multiplies by X^w, and the digit mask keeps exactly the
required low degrees. Part two recomputes every entry of
the two exact arrays in admitted.json.

## Competitors changing every field

Write q=2^50 and Q=q^g. For nonzero e, the cube
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
directions e_i in all g fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the g translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the K
parity bits can be evaluated by expanding the
g affine indicators. Multiplication by their
common point mass 2^g/Q leaves 1/Q times a sum
of 2^g terms. The empty mask term is the
parity ball volume. Each nonempty term has an
absolute bound equal to the absolute signed sum
of the parity ball against the combined mask of
its selected parity components.

For a K-bit mask of weight w>0 the signed ball
sum is, up to translation sign, the coefficient
K_l(w) of X^l in (1-X)^(w-1)*(1+X)^(K-w).
Indeed the polynomial of shell signs has the same
form with exponent w on (1-X), and summing the
prefix divides by 1-X. Let N=K-1 and
gamma=N-2*(w-1). Differentiating this polynomial F gives
(1-X^2)F'=(gamma-NX)F. The coefficient recurrence is

    (l+1)*K_(l+1) = gamma*K_l-(N-l+1)*K_(l-1),

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
H_k(l,c) to be the largest sum of c absolute
magnitudes in this list, sorted by |K_l(w)| with
the displayed multiplicities. The computation takes
full groups of largest magnitude then a partial
group, giving the maximum over any at most c
directions. The available multiplicities exceed the
capacities used.

Fix the true noise. In a chosen k-mask term the
relevant k systematic fields are the predecessors
of the selected parity fields. For a combined
systematic error weight a on them, there are at
most binom(50*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(50*(g-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(50*(g-k),b) *
    H_k(t-a-b, binom(50*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(g,k) selections of masks, and add
the n-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_t. The expected
number of full-change competitors within total
t for each fixed true noise is at most D_t/Q.
The checker verifies D_t is nondecreasing through t=T.

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
over messages. Convolve this by the counts of outside noises of weight t-w.
This gives the upper count at each exact global noise weight.
Sum over the g cyclic starts and all proper run lengths.
The remaining calculations use exact integer or rational
arithmetic.

## Occupancy after exact-weight removals

At each total true noise weight v upper-bound the mean count
of pairs removed per true message. Remove the pair if its
true anchor or child path fails admission or if a proper-run
competitor is closer, or is tied and lexicographically earlier.
The fallback may still supply extra candidates.

For anchor systematic weight d, anchor parity weight w and
child weight z=v-d-w, their disjoint counts are multiplied.
The anchor factors are binom(50,d) and binom(50,w).
For d>R every child choice of weight z is excluded. Otherwise
the number excluded is binom(n-100,z)-A_(d,z).
Sum these products over nonnegative d,w,z with d+w+z=v.
This is the exact count of admission failures at each weight.
Add the proper-run upper count described above, round upward
to an integer, and clip at binom(n,v). Denote this upper
mean count by U_v. Repeated removals can only enlarge it.

At t<=T let V_t be the length-n Hamming ball volume and let
H_t=V_t-sum_(v<=t) U_v. At least Q*H_t ordered pairs of
true and received words remain, where Q=2^K. Since any two
remaining distinct words at the same receiver differ
in all fields, their ordered pair count is at most V_t*D_t.
This follows by summing the uniform full-change bound over
all true noise at distance at most t. All remaining words
are in the actual decoder list.

For the count m_y of remaining words at receiver y and any
integer k>=1,

    k*(k+1) * 1[m_y>0] >= 2*k*m_y-m_y*(m_y-1).

The difference at positive m_y is (m_y-k)*(m_y-k-1),
and it vanishes at zero. The number F_t of received words
where the decoder finds a distance at most t thus obeys

    F_t/Q >= (2*k*H_t - V_t*D_t/Q)/(k*(k+1)).

For H_t>0 use k=1+floor(V_t*D_t/(Q*H_t)). Floor the
rational lower bound, clip at zero, and take its running
maximum G_t. Use the zero bound when H_t=0. Monotonicity of
F_t justifies the running maximum.

At exact decoded distance t, the success probability for
a particular received word and its selected message is
9^(n-t)/(Q*10^n). The decoder keeps a least-distance
candidate in its list. Summation by parts with the
decreasing weights, and ignoring distances beyond T,
gives success at least

    sum_(t=0..T) (G_t-G_(t-1))*9^(n-t)/10^n,

with G_(-1)=0. The calculation proves an exact conservative
component bound directly from the verified admission counts.

## Combining the components

Write b3 and b2 for the exact upper error bounds from the
component calculations with g=3 and g=2, respectively.
Within the uniform 1000-bit message their payloads are
independent, and all their channel noises are independent.
Each reported component value is a deterministic function
of its own observations. Every local success event is
therefore independent of the others. The block decoder
succeeds when all seven component decoders succeed, giving
the exact upper bound 1-(1-b3)^6*(1-b2). Compute this rational
before rounding upward to the declared endpoint.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
For each component its 2g received fields are copied to padded little-endian
one-word fields in arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
Words SAVED_B through SAVED_B+19 retain all twenty decoded field values.
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
- Each component configuration has its own budget_base, one word per anchor depth 0 through R holding C(d).
- The configuration profile24 and profile1 arrays have one word per orbit (width 24) or syndrome
  (width 1) giving its complete-profile type number.
- Its cap_offsets[i-1] is the base for emitted cell i. Store
  E_(i,p)(r) at offset p*851+r for every type p and state r.

The files construction.py, qr.py and dynamics.py specify every fixed datum. The exclusive end of RAM is the value MEMORY_END in
construction.py, within the model's 2^32-word limit. All table addresses are distinct and the
program's temporary words precede these tables.
The two configurations share the field, QR, and repetition banks.
The two disjoint configuration-specific regions follow the QR
row bank. Their literal budget, profile and cap addresses are
specified by component_layout.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All J DFS frames use RAM, and the cell and
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
For component c, let f be the number of fields in the preceding
components. Extract its g systematic fields starting at bit 50*f
and its g parity fields starting at bit 1000+50*f. Store each
extracted field in the corresponding padded word of Y or P.
Load, popcount and add each of those words for the distance of
the all-zero incumbent, then clear the g incumbent B words.
Copy Y0 to X0 in two primitives, set Depth to zero and
store NEXT[0]=0. This initial setup uses at most 100 primitives
per component, since g<=3. Each cell and anchor frame writes
all its slots before reading them.

When packing a field from short pieces, clear an accumulator
by a free move, then insert each piece by a left shift and OR
at its literal bit offset, and store the result.

Emit the entire component traversal seven times with the
literal configuration chosen at emission. After each return
copy the g padded incumbent words to SAVED_B+f in 2g operations.
Then emit the sixteen output words in order. Gather each from
at most three successive padded source words of width 50.
Extract the required pieces and insert them into the cleared
output accumulator at their literal offsets, and store it.
Each word costs at most 3*(6+2)+1. The final word contains
only forty data bits and its upper twenty-four bits stay zero.
All outer initialization, component-result copies, and output
gathering costs at most 7*100+40+16*25=1140. The calculation
uses a conservative outer allowance of 50,000.

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
word, XOR and store to AP. Set Remaining by adding Depth to the literal configuration budget_base and loading
that word. Move Distance=Depth and jump to the first
CELL_INIT below.
This uses one cube and at most 250 further primitives:
the four cells use at most 4*(6+2+4), and finishing AP
uses 3, two budget operations, and one jump.

### Cell traversal

Emit a literal block for each of the J cells in field and within-field order.
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
        idx := add(profile1,s)
        p := load(idx)
        j := mul(p,851)
        j := add(j,Remaining)
        idx := add(cap_offsets[DepthCell-1],j)
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
        needed := add(cc,J-DepthCell)
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
        run FIELD after each of the nonfinal child fields
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
    idx := add(profile24,orbit)
    p := load(idx)
    j := mul(p,851)
    j := add(j,Remaining)
    idx := add(cap_offsets[DepthCell-1],j)
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
the anchor parity to Distance by loading AP and C_(g-1) once, XORing, popcounting and adding.
The residual AP+C_(g-1) is P0+L(X0)+C_(g-1). Distance already includes
the anchor systematic distance and all g-1 visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all g systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 200 primitives.
## Exact worst-case operation bound

Number the J literal cells by i=1..J. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 67 at width 24 and 47 otherwise, and g_i be 80 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
T_i is CubeOps+180 at every multiple of four before J,
and CubeOps+200 at i=J. It is zero otherwise.

For all possible states r=0..850 define in decreasing cell order

    F_(J+1)(r)=0
    U_(i,p)(s)=sum_(c<=s-(J-i))
                 A_(p,c) * (g_i + T_i + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,850]: U_(i,p)(s)<=F_i(r)-h_i}.

The nonnegative recurrence is nondecreasing in r by induction
and U is nondecreasing in s. The maximum defining E exists
because r itself is allowed. Its search returns the rightmost
valid index, including equality. These are the fixed CAP table entries.
For any actual state r and profile p, the accepted attempts
are exactly the rows up through their recorded cost threshold
E_(i,p)(r)-(J-i). Their continuation states E_(i,p)(r)-c
belong to [0,850] because costs are positive and successful
attempts satisfy the reserve. At a successful attempt the
remaining program cost is bounded by g_i+T_i+F_(i+1) at that state.
Consequently the actual entire cell traversal costs at most
h_i+U_(i,p)(E_(i,p)(r))<=F_i(r). This proves the bound by
downward induction without distributional hypotheses.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 250+CubeOps. The per-component upper count at its literal configuration is

    W_g = 1000 + sum_(d=0..R) binom(50,d) *
                  (44 + 250 + CubeOps + F_1(C(d))).

The added 1000 is harmless padding of each local bound.
The mixed decoder has the conservative all-input upper count
W=50000+6*W_3+W_2=8361251488734102, including the separately bounded outer
work. The declared exact recurrence checks W<=2^(52893/1000).

## Declared computations

Part zero checks the shared QR spectra, orbit and field
identities, and literal layout bounds. Part one recomputes
the exact adaptive cap and primitive-count recurrences for
both component configurations. Part two recomputes every
stored admitted child coefficient using those caps and
the verified profiles. Part three combines exact noise-shell
removals, cumulative occupancy and the independent product
formula above. It verifies B=13699779932857317/15625000000000000 and the exact rational
work exponent from the same fixed decoder. No physical
materialization of the largest table is required.
