# Five independent four-field cyclic codes with adaptive decoding

## Encoder and decoder

Use the polynomial basis of F2[X]/(P), where
P=X^50+X^9+X^2+X+1. The declared field check applies
the Rabin criterion with prime divisors 2 and 5 of 50.
Split the 1000-bit message into five 200-bit parts, each split
into four 50-bit words b0 through b3 in increasing bit order.
Within each part define p_i=L(b_i)+b_(i-1)^3 with cyclic
indices modulo four. Transmit the five parts' systematic words
in order, then their parity words in the same part and field order.
Each part is a 200-bit systematic encoder with 400 channel
coordinates. The 2000-bit encoder is injective by its systematic
prefix. The decoder below processes each component on its own
400 received coordinates and concatenates the five decisions.
Reused scratch and tables do not introduce dependence between
decoded component values.

In each field L is block diagonal with two 24-bit cells and two
width-1 repetition cells in increasing bit order. For a width-24
input u, the low 23 parity bits are the binary polynomial
remainder of u*X^23 modulo the degree-23 G=9205487 (in
low-coefficient-first order). The high parity bit is the parity
of the sum of Hamming weights of u and the low parity bits.
The width-1 parity map is the identity. These maps are linear. Construction.py fixes the integer layout.

The decoder enumerates all anchors x0 of distance d<=11
from received Y0. For each use entry C(d) at index d of construction.py BUDGETS as the initial
remaining-cost state. Visit the child fields j=1 through 3 in order. In each field,
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

Depth-first search the 12 child cells, retaining a prefix of length i
only at cost at most E_(i,p)(r)-(12-i), where r is
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

## Exact child-cost distributions

On the true path the predecessor cube equals the true cube and
subtracting it leaves the original independent observations.
For the channel with crossover probability 1/q, where q is one
of 10 and 12, a particular cell noise of weight w in width m
has probability (q-1)^(2m-w)/q^(2m). The decoder always uses
the nominal cost defined from the weight base 9. The auxiliary
channels only bound admission probabilities.

For each width m, group syndromes by their complete noise profile and
number the distinct profiles in lexicographic order. A type p
with multiplicity M_p supplies k_(p,w) retained patterns of true weight w
in each syndrome, where retention takes the least weights until
the row cap is exhausted. Let A_(p,c) count its retained entries of nominal
cost c. The following local numerator counts true patterns across
all syndromes of the type under the 1/q channel:

    H_(q,p,c) = M_p * sum_(w of cost c) k_(p,w) * (q-1)^(2m-w).

Orbit permutations preserve weights, so a cap that cuts a distance tie
has exactly these retained counts for every syndrome in the orbit.
The RAM maps an orbit index to its representative's type
and a small-cell syndrome directly to its type.
The program data generator defines each such entry by enumerating
the complete syndrome noise profile. The tested histogram and
the mapping use the same profiles and lexicographic numbering.

For cell i and nominal remaining state r, the program looks up
the effective cap E_(i,p)(r) of its actual type p. The operation proof
defines these caps on all states 0..850 by the exact recurrence below.
With s=E_(i,p)(r), the entry of cost c is accepted precisely when
c<=s-(12-i) and the entry is in the row. The next state is s-c.
This description also applies to the true path: its syndrome
at each cell depends only on that cell's true noise, since on that
path the cube adjustment is correct. True cell noises are
independent under each channel, including conditional on previous
true noise cells and the anchor. Thus for each q the exact
conditional suffix-admission probabilities satisfy

    P_(13,q)(r)=1
    P_(i,q)(r)=q^(-2m) * sum_p sum_(c<=E_(i,p)(r)-(12-i))
                        H_(q,p,c) * P_(i+1,q)(E_(i,p)(r)-c).

Compute lower integers with scale S=10^60 using precisely this
recurrence and floor each division by q^(2m). The starting integer
for the empty suffix is S. All factors are nonnegative, so
the rounded integers are at most S times the probabilities by
downward induction. For anchor depth d, the complement through one
is an upper bound on the probability E_d that the true child path is
absent, when evaluating state r=C(d) at the first cell.
The arrays in tails.json contain this upper bound in units 1/S
for q=10 and q=12 at all depths 0..11. Declared parts 1 and 2
check every entry. This event depends only on child noise and
the stated depth, which fixes its initial state.

## Competitors changing all four fields

Write q=2^50 and Q=q^4. For nonzero e, the cube
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
directions e_i in all four fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the four translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 200
parity bits can be evaluated by expanding the
four affine indicators. Multiplication by their
common point mass 16/Q leaves 1/Q times a sum
of 16 terms. The empty mask term is the
parity ball volume. Each nonempty term has an
absolute bound equal to the absolute signed sum
of the parity ball against the combined mask of
its selected parity components.

For a 200-bit mask of weight w>0 the signed ball
sum is, up to translation sign, the coefficient
K_l(w) of X^l in (1-X)^(w-1)*(1+X)^(200-w).
Indeed the polynomial of shell signs has the same
form with exponent w on (1-X), and summing the
prefix divides by 1-X. Let N=199 and
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
binom(50*(4-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(50*(4-k),b) *
    H_k(t-a-b, binom(50*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(4,k) selections of masks, and add
the 400-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_t. The expected
number of full-change competitors within total
t for each fixed true noise is at most D_t/Q.
The checker verifies D_t is nondecreasing through t=59.

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
over messages. The computation multiplies this by the
counts of outside noises of weight t-w and by the
channel probability 9^(400-t)/10^400.
It sums the resulting error mass with alpha_t. All
operations in the conditional count and its application
are exact integer or rational arithmetic.

## Occupancy of the actual candidate list

At distance cutoff t<=59 consider all pairs of true word and
received word with true distance at most t. Remove any pair
for which the true word fails the specified anchor or child
admission tests, and any pair for which a single-run competitor
is closer or is tied and lexicographically earlier. The fallback is allowed even when it
would have been removed. Let L_t count removed pairs, Q=2^200,
and V_t be the 400-bit Hamming ball volume through t. Then
Q*V_t-L_t pairs remain. Two distinct retained candidates at the
same received word must differ in all four fields, since otherwise
the farther would be removed, or the lexicographically later on a tie.
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
success is 9^(400-t)/(Q*10^400). Thus the success through the
cutoff is obtained by summation by parts with w_t=9^(400-t)
and w_60=0, by multiplying F_t/Q by (w_t-w_(t+1))/10^400
and summing. Ignoring success outside the cutoff only decreases
the lower bound. A removed pair with true noise weight v can
be charged only at t with i(t)>=v, and alpha_(i(t))<=alpha_v
because D_t is nondecreasing. Its total charge is at most
alpha_v*w_v/(Q*10^400) by telescoping.

The error bound therefore consists of the unremoved baseline

    1 - sum_(t<=59) (U_t-U_(t-1))*w_t/10^400,

plus the proper-run terms for all four starts, plus the admission
failure charge weighted by alpha at the true total noise weight.
Here U_(-1)=0. For each run length a=1..3 use J_(a,v) above, multiplied
by alpha_v and summed through 59, and charge all four starts. An upper
bound on the weighted admission charge may count both removal
reasons even when they occur on the same pair.

## Weight-conditioned admission charge

Write d for the true anchor systematic weight on 50 bits.
Write w for its independent parity noise on the other 50 anchor
bits and Z for the child total noise on 300 bits. The total weight
is d+w+Z, with the three noise blocks mutually independent.
If d>11 the anchor is absent whatever the other noises; its
exact weighted charge is

    sum_(d>11) Pr[Bin(50,1/10)=d]
        * sum_(z<=59-d) Pr[Bin(350,1/10)=z] * alpha_(d+z).

Use only d<=59 because the weighted charge vanishes otherwise.

If d<=11, failed admission is exactly E_d under the child's
nominal noise. For q=10,12 the likelihood ratio of a child
noise pattern with weight z under 1/10 to its probability under
1/q is

    R_q(z)=q^300*9^(300-z)/(10^300*(q-1)^(300-z)).

It is nondecreasing in z, since its step multiplier is (q-1)/9.
Hence for arbitrary E_d,

    Pr_(1/10)[E_d and Z<=z] <= min_q R_q(z)*Pr_(1/q)[E_d].

This inequality covers both failure to find a required QR row entry
and failure of the cost test. It does not assume a monotone relation
between noise weight and normalized cost.

Put alpha_60=0 and define

    G(s)=sum_(w<=59-s) Pr[Bin(50,1/10)=w]
                            *(alpha_(s+w)-alpha_(s+w+1)).

For fixed d,w, summation by parts over z turns the truncated
weighted child charge into nonnegative multiples of the cumulative
probabilities Pr[E_d and Z<=z]. Independence of w permits summing
first over it, producing G(d+z). The child contribution is thus at
most the sum for d<=11 and z<=59-d of

    Pr[Bin(50,1/10)=d] * G(d+z) * min_q R_q(z)*Pr_(1/q)[E_d].

The verifier computes the exact anchor term. For the child term
it rounds G up to denominator 10^60, multiplies the verified tail
ceilings by the exact ratios, and rounds these products up to
denominator 10^60 before minimizing. Each replacement is upward,
and both weights are nonnegative. It then sums the resulting
integer numerators against the exact anchor binomial masses.
Part 3 checks this weighted error bound together with the full
change and proper-run contributions.

## Combining the components

Let b be the exact upper bound given by the four-field occupancy
and miss calculation above for one component. Component inputs
are independent uniform 200-bit messages. The respective sets
of 400 channel noises are independent. Each component result
is a deterministic function of its own observations. Each success
event therefore has probability at least 1-b and the five events
are independent. The final decoder succeeds exactly when all five
components succeed, so its block error is at most 1-(1-b)^5.
The exact arithmetic computes that quantity before rounding up.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
For each component the eight received fields are copied to padded little-endian
one-word fields in arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, six frame words per child cell, and cube scratch.
Words SAVED_B through SAVED_B+15 retain the first four component results.
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
  E_(i,p)(r) at offset p*(851)+r for every type p and state r.

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
can be discarded. All 12 DFS frames use RAM, and the cell and
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
For component k=0 through 4 emit a literal copy from
systematic starting bit 200*k and parity starting bit 1000+200*k.
Each group contains four fields of 50 bits. Use the extraction
routine for eight padded words of Y and P, then store each.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 24. Store zero in the four
padded incumbent words. Copy Y0 to X0 in two primitives,
set Depth to zero and store NEXT[0]=0. This setup costs at most
87 primitives per component. Each visited cell and anchor frame
writes every slot it later reads.

For packing a field from known short pieces, clear one
register by a free move and insert each piece with a left shift
and OR at its literal bit offset. Store the accumulator at
the field's word address.

Emit the traversal below five times with the corresponding
literal input offsets and fresh setup. After each of the first
four returns, copy the four padded incumbent words B to
SAVED_B+4*k using eight primitives. After the last return,
the sixteen output words gather the first 800 bits from SAVED_B
and the last 200 bits from the current B.
Each output word receives pieces from at most three successive
padded source words, since the source fields have width 50.
Extract each literal piece and insert into a zeroed output
accumulator; this costs at most 3*(6+2)+1 including the output
store. The last output word has its high 24 bits zero. Combined
non-traversal work including setup, transfers and stores
is at most 5*87+4*8+16*25 <1000.

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

Emit a literal block for each of the 12 cells in field and within-field order.
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
        j := mul(p,851)
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
        needed := add(cc,12-DepthCell)
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
        run FIELD after each of the two nonfinal fields
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
    j := mul(p,851)
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
the anchor parity to Distance by loading AP and C3 once, XORing, popcounting and adding.
The residual AP+C3 is P0+L(X0)+C3. Distance already includes
the anchor systematic distance and all three visited pairs,
so this is the exact full distance. Compare with BestD, and
on strict improvement store all four systematic X words
in B and update BestD. Including its packing, the noncube
FINAL charge is at most 200 primitives.
## Exact worst-case operation bound

Number the 12 literal cells by i=1..12. For the cell's width,
let A_(p,c) count its retained row entries of cost c for each
profile type p. This is the count used in the true-path proof.
Let h_i be 67 at width 24 and 47 otherwise, and g_i be 80 at
width 24 and 28 otherwise. The setup bounds include the
type lookup, effective cap lookup, and one stopping attempt.
T_i is CubeOps+180 at i=4,8 and CubeOps+200 at i=12.
It is zero at all other indices.

For all possible states r=0..850 define in decreasing cell order

    F_13(r)=0
    U_(i,p)(s)=sum_(c<=s-(12-i))
                 A_(p,c) * (g_i + T_i + F_(i+1)(s-c))
    F_i(r)=h_i+max_p U_(i,p)(r)
    E_(i,p)(r)=max {s in [0,850]: U_(i,p)(s)<=F_i(r)-h_i}.

The nonnegative recurrence is nondecreasing in r by induction
and U is nondecreasing in s. The maximum defining E exists
because r itself is allowed. Its search returns the rightmost
valid index, including equality. These are the fixed CAP table entries.
For any actual state r and profile p, the accepted attempts
are exactly the rows up through their recorded cost threshold
E_(i,p)(r)-(12-i). Their continuation states E_(i,p)(r)-c
belong to [0,850] because costs are positive and successful
attempts satisfy the reserve. At a successful attempt the
remaining program cost is bounded by g_i+T_i+F_(i+1) at that state.
Consequently the actual entire cell traversal costs at most
h_i+U_(i,p)(E_(i,p)(r))<=F_i(r). This proves the bound by
downward induction without distributional hypotheses.

The anchor subset traversal uses at most 44 primitives
per visit as counted above. Parity residual setup and the
anchor cube cost at most 250+CubeOps. Hence a worst-case
all-input integer bound for all five traversals and their layout is

    W = 1000 + 5 * sum_(d=0..11) binom(50,d) *
                 (44 + 250 + CubeOps + F_1(C(d))).

The exact integer recurrence in the declared verifier
computation checks W against the declared rational exponent.

The declared computation gives the following ceilings with the
common denominator 10^18 for the probability components.

    Unremoved baseline error <= 214118756748517276/10^18
    Single-start proper-run charges <= (189957540195583, 455973364326611, 93057657901850)/10^18
    Anchor admission charge <= 1176167154322930/10^18
    Child admission charge <= 38612321118644217/10^18
    B = 773356200462581483/1000000000000000000
    W = 2520497841559917180
    W^(1000) <= 2^(61129).
    Exclusive RAM end = 3686885648.

The exact sum charges each of the 4 starts in a component.
The component error is at most 256863199271180589/10^18.
The computation applies the independent product bound to the exact
component value, then rounds the combined result upward. Separately rounded
component displays may sum slightly above the combined ceiling.
