# Cyclic cube parities with an integer inner budget

## Code and decoder

Use five 200-bit information fields b_0 through b_4
in little-endian order, with coordinates in the
polynomial basis of K=F2[X]/(P), where

    P = X^200 + X^30 + X^19 + X^11 + 1.

The declared computation verifies X^(2^200)=X modulo
P and gcd one between P and X^(2^k)+X for k=40,100.
An irreducible factor of degree d would have d divide
200 because its generating root is fixed by the 200th
Frobenius power. Any proper divisor of 200 divides
one of 40 and 100. Thus the two gcds prove irreducibility.

Transmit all five b_i followed by their five parity
fields in the same index order, with

    a_i = b_i xor b_(i-1 mod 5)^3.

This is an injection carrying exactly 1000 bits in
2000 binary channel uses. All b_i are independent and
uniform under uniform information.

The decoder begins with the all-zero candidate.
Write Y_i and P_i for the received systematic and
parity fields. Choose each i in turn as an anchor,
in order 0,1,2,3,4. Traverse all subsets of up to
R=24 flips to Y_i. For anchor weight u the
following inclusive bands determine the total
remaining integer budget K(u), in bits.

    0 <= u <= 13: K(u) = 266
    14 <= u <= 16: K(u) = 258
    17 <= u <= 19: K(u) = 248
    20 <= u <= 24: K(u) = 233

For each anchor candidate visit its other four
fields in cyclic order. At a child field j the
predecessor candidate x_(j-1) is known. Compare
Y_j with P_j xor x_(j-1)^3. Suppose they disagree
at d positions. Any proposed child x_j is specified
by flips relative to Y_j. If h flips are on agreeing
positions, the child costs

    5+d+6*h

units. Keep prefixes whose total cost is at most the
anchor budget minus five per not-yet-chosen field.
For each child this permits all disagreement flips
and all agreement subsets through a specific cap.
The reserve ensures possible remaining baseline
costs fit, so it does not exclude any complete path
within budget. A path of four children within the
total budget is scored by its full codeword distance.
Return the earliest minimum over both the fallback
and the five orientation searches. The precise order
is the subset tree order in the scalar program.

## Sufficient admitted shells

Conditional on a true predecessor, the two received
observations of a child differ from its true field
by the noises in its own systematic and parity fields.
If those noises have d disagreements and h double
flips, their group contributes weight d+2h. Given
systematic weight u, the corresponding number of
patterns is

    binom(200,d)*binom(200-d,h)*binom(d,u-h).

Along the true path a disagreement has no constraint,
and h is the number of needed agreement flips.
Thus its four children have exact total cost
20 plus the sum of their values d+6*h. The
reserve in the decoder permits those children if
this total is within the anchor's budget.

For the lower count choose a preferred true anchor
for each noise pattern as follows. Find the first
band occupied by any of the five systematic noise
weights, and take the first index in that band.
If no band is occupied, this pattern is omitted.
Count the pattern only if this particular preferred
anchor's true path satisfies the budget. The real
decoder searches all anchors, so this suffices for
admission. The selection here depends only on the
noise pattern, and is used for the certificate.

Fix a band with endpoints low,high and budget K.
For the preferred anchor to belong to it, its own
systematic noise weight is in that band. Each earlier
index has systematic noise weight greater than high.
Each later index has systematic noise weight at
least low. These are exactly the conditions for the
specified band and index, and all five possible
preferred indices give disjoint events. The earlier
and later groups are the same four children of that
anchor, in some cyclic order. For actual true paths
their total cost is independent of the traversal order.

Make two polynomials of one child group in variables
(C,H). In each polynomial the monomial for d
disagreements and h double flips is C^(d+6*h)*H^h.
Its coefficient is the preceding binomial expression
summed over systematic weights at least low for
polynomial A, and greater than high for polynomial B.
For the five possible preferred-index placements,
the combined polynomial of the other four groups is

    A^4+A^3*B+A^2*B^2+A*B^3+B^4 = (A+B)*(A^3+B^3)+A^2*B^2.

The budget truncates its C degree to K-20. For a
term C^c*H^h, d=c-6*h and its total noise weight
is d+2*h. Combine these child counts by ordinary
weight convolution with the preferred anchor's own
200 systematic and 200 unrestricted parity bits.
The resulting 231 coefficients count exactly the
sufficiently admitted patterns of total weight
0..230 for that preferred band. The four bands'
sufficient events are disjoint.

In the experiment use stride s=1+floor((K-20)/6)
and encode C^c*H^h at degree s*c+h as a packed
univariate coefficient. Any term with h>=s has
c>=6h>K-20, so cannot alias a retained degree.
Terms with c>K-20 lie beyond all retained degrees.
For each term with c<=K-20, its true noise weight
is at most c. Therefore every retained coefficient
in a product of at most four group factors or in
the displayed sum is at most
5*binom(1600,K-20). This follows by
dominating by all unrestricted noise patterns of
that weight and using monotonicity of the binomial
coefficient in the indicated prefix. The declared
ranges have K-20<=800. The checker verifies
this bound is below the coefficient radix 2^1120.
Integer multiplication therefore computes these
truncated polynomials exactly.

The four declared runs each reproduce 231 coefficients
for one distinct band and compare them to the corresponding
array in shells.json. The final run sums all four checked
arrays and evaluates the full exact certificate. Define E_t to be these sufficiently
admitted noise patterns through weight t, and J_t
their cumulative number. The subsequent argument
only needs E_t to imply actual decoder admission.

## Competitors changing all five fields

Write q=2^200 and Q=q^5. For nonzero e, the cube
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
directions e_i in all five fields. Its parity
difference in field i is e_i plus the cube
difference in field i-1. Thus, for fixed directions,
the five translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 1000
parity bits can be evaluated by expanding the
five affine indicators. Multiplication by their
common point mass 32/Q leaves 1/Q times a sum
of thirty-two terms. The empty mask term is the
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

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(200*h,w).

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
most binom(200*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(200*(5-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(200*(5-k),b) *
    H_k(t-a-b, binom(200*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(5,k) selections of masks, and add
the 2000-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_t. The expected
number of full-change competitors within total
t for each fixed true noise is at most D_t/Q.
This estimate may count candidates outside the
actual search. The checker verifies D_t is
nondecreasing through t=230.

## Conflicts from proper runs of changes

For removal only, order codeword candidates by increasing
(distance to y, numerical message), where the message is
the 1000-bit little-endian integer. This is a strict total
order independent of the decoder's own ties. A nonempty
proper subset of changed systematic fields on the cycle
splits into maximal consecutive runs, each of length in
{1,2,3,4}. A run changes its own systematic and parity
fields and the next parity after it. These affected
coordinates are disjoint between runs.

If a proper-change word is better under this order, either
some run strictly improves the distance, or all run
distance changes are zero and the word's numerical message
is smaller. For the latter case the run containing the
most significant differing message bit is itself better.
If there is a negative distance change, its run suffices
regardless of the numerical order. Otherwise equal total
distance requires every local distance change to be zero.
Thus some single-run competitor is better in the stated
order. It suffices to remove a candidate when this occurs.

For one fixed run start and length a, the first
systematic field's nonzero direction is repeated
in its parity, since its predecessor does not
change. Each of the other a parity differences is
an affine trace distribution from a distinct true
field, with a deterministic systematic offset
depending only on directions. Conditional on a
fixed direction vector, these are independent and
have joint point mass at most 2^a/q^a. For a fixed
true noise on the affected coordinates, a possible
competitor error on all systematic fields fixes
all directions. The repeated pair must preserve
its disagreement set. If this set has d positions,
the weight polynomial for all possible errors
on the affected coordinates with this same set is

    P_(a,d)(X) = 2^d*X^d*(1+X^2)^(200-d)*(1+X)^((2a-1)*200).

The repeated pair gives the first three factors.
The other (2a-1)*200 positions comprise the
remaining systematic and varying parity fields.
This polynomial counts all candidate patterns on
the pair, even ones implying a zero direction,
so is an upper envelope. For a fixed disagreement count d and true local
noise weight w, write V_(a,d)(w) for the cumulative
coefficients through w. Conditional on any fixed
true noise in that stratum, the expected number of
run competitors at local distance at most w is
at most 2^a*V_(a,d)(w)/q^a, by the parity point
mass bound. For strictly closer the same bound
uses V_(a,d)(w-1), with the count zero at w=-1.

For equally distant ordered pairs, average over
uniform true messages and noise in this stratum.
Swapping the true message and its tied competitor
preserves channel likelihood, the disagreement set
of the repeated pair, the unchanged outside noise
and the local distance w. Their codewords differ
on the same single run. Exactly one direction of
the swap is better in numerical message at this
tie. Thus half the expected tie count is charged.
The expected number of preferred run competitors
in this stratum is at most the average of the two
cumulative bounds just proved. By Markov and
clipping, the preferred-run conflict probability
has an upper numerator

    min(q^a, 2^(a-1)*(V_(a,d)(w)+V_(a,d)(w-1)))

over denominator q^a. This use of clipping is on
the uniform conditional average within the stratum.
The other noisy bits outside the affected coordinates
are independent and do not affect the comparison.

For each w construct S_(a,w) by multiplying this
upper numerator by binom(200,d)*[X^w]P_(a,d) and
summing over d. The total number of full noise
patterns weighted by this numerator at noise
weight v is bounded by

    T_(a,v) = sum_(w=0..v) S_(a,w)*binom(2000-(2a+1)*200,v-w).

For a fixed run start and length a the probability
of a preferred run conflict at weight v is at most

    T_(a,v)*9^(2000-v)/(q^a*10^2000).

The verifier computes the unshifted coefficients
(1+X^2)^200*(1+X)^((2a-1)*200) and successively
divides by 1+X^2 as d increases, using the recurrence
new[i]=old[i]-new[i-2]. It uses the exact polynomial
prefixes through weight 230. There are
five possible run starts, so the union bound over all
lengths uses five times the sum of these bounds.

## Occupancy after removal

Among the admitted pairs (u,y) with noise in E_t,
remove any candidate for which there exists a
run competitor better in the stated order. Let L_t be the
cumulative number removed. If two distinct
retained candidates at a common y differed on a
proper subset of fields, take the one worse in
the removal order. The other is a better proper
competitor, so a single run would force removal.
This is a contradiction. Thus all
ordered distinct pairs of retained candidates at a
common y are full changes. The double sum of
their numbers over y is at most J_t*D_t by summing
the full-change conditional upper bound over the Q
true words at each fixed admitted noise.

Let F_t count received words having any actual
decoder candidate of distance at most t, and m_y
count retained candidates. Each retained candidate
is in the actual search. For integer k>=1,
k(k+1)1[m_y>0] >= 2k*m_y-m_y*(m_y-1), since
(m_y-k)*(m_y-k-1)>=0 for positive m_y and both
sides vanish at zero. Applying sum_y m_y=Q*J_t-L_t
and the ordered pair upper bound gives

    F_t/Q >= G_t - alpha_t*L_t/Q
    G_t = floor(J_t*(2*k_t*Q-D_t)/(Q*k_t*(k_t+1)))
    k_t = 1+floor(D_t/Q),  alpha_t=2/(k_t+1).

Because D_t is nondecreasing, alpha_t is
nonincreasing. Put U_t=max_(i<=t) G_i and choose
an index i(t) attaining this maximum. Monotonic
F_t gives F_t/Q>=U_t-alpha_(i(t))*L_(i(t))/Q.

For each received word, the decoder chooses one
candidate of smallest distance. Its likelihood,
divided by Q, equals that word's joint success
contribution. Let w_t=9^(2000-t) for t<=230 and
w_231=0. The joint success through the cutoff
is obtained by summation by parts, which uses
the nonnegative differences w_t-w_(t+1) on the
cumulative F_t/Q and denominator 10^2000.
Substitute the lower bound above. A removed
pair of true weight v contributes to its loss
only when i(t)>=v, which implies t>=v and
alpha_(i(t))<=alpha_v. The total coefficient over
those t is at most alpha_v*w_v by telescoping.
After normalizing by Q*10^2000 this charges at most
alpha_v times the probability of a removed pair at
weight v. The run bounds therefore give overall
error at most

    1 - sum_(t<=230)(U_t-U_(t-1))*w_t/10^2000
      + 5 * sum_(a=1..4) sum_(t<=230) alpha_t*T_(a,t)*w_t/(q^a*10^2000),

where U_(-1)=0. Omitted true-success terms are
nonnegative. The declared computation checks this
quantity is at most 11123976282903969/125000000000000000 using exact arithmetic.

## RAM representation and field arithmetic

Each field is a padded little-endian array of k=4
words with h=8 live bits in the last word.
Construction.py specifies every address. Starting
at word 48 put the six groups of five field arrays
Y, P, X, C, I, B, in that order, then SQ. Y and P
store the received systematic and parity fields.
X and C hold the current candidates and their
cubes. I holds the cubes of the corresponding
received systematic fields. B holds the incumbent
message, initially zero. Each of nine trees has
a disjoint NEXT and LAST array of 201 words. After
those put the two 200-word coordinate lists for
each of the four child levels, then the two
201-by-k saved-cube arrays for each of those levels.
The multiplication scratch holds sixteen k-word
rows. Then come the fixed reduction, square, depth budget and affine tables. All are consecutive with
exclusive end 2516619929, within the 32-bit address
limit. The specified addresses do not alias the
input or output prefix.

All fixed ranges and choices in the following
program expand at code generation. Macros inline
their primitive lists. Literal addresses, constants
and register moves are free. Dynamic tests, jumps,
loads, stores and arithmetic are charged as shown.
The data generators fix all table words before
decoding. Their polynomial squares interleave bits
with zeros and fold high terms by the sparse modulus.

For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^200. This fits in one word since the
largest possible degree is 33. To multiply a field
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
c*X^200 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=25 byte slots, emit a
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
Move R registers to zero. For each of n0=50 nibble
slots of the source, in decreasing order, load its
word, shift right, AND with 15, multiply by k and
add the multiplication base. Apply SHIFT_4 followed
by adding the indicated row to R with one address
add, load and XOR per word. Horner's rule yields
the field product. Store R to the designated cube
array. The multiplication cost is

    Mmul=3k+3*(3k+3+k)+44k+n0*(5+3k+3+3k)+k.

Both parts together have cost C0=S+Mmul and give
the exact field cube with zero high padding. The
high padding of the source supplies zeros in the
last partial byte or nibble.

## Affine updates and subset tree

At a non-anchor tree descent flipping coordinate i,
save the k current cube words in row depth of the
tree's saved-cube array. Use mul and add for the row
base and an add, load and store per word. If x is
the parent field and e=X^i, its cube change is

    (x+e)^3+x^3=e*x^2+e^2*x+e^3.

For direction i the fixed affine slices cover 12
successive chunks of 18 bits, up to position 215,
with positions 200 through 215 set to zero. The entry
d of slice l XORs the linear columns
e*X^(2j)+e^2*X^j for all its selected live bits j.
It includes e^3 exactly for slice zero. Entry
word h0 is at

    AFFINE_BASE+i*AFFINE_STRIDE+l*262144*k+d*k+h0,

with AFFINE_STRIDE=12*262144*k. The generator
computes this entry at every fixed address.

Compute the direction base with mul and add, and
load the current cube as k accumulators. For
each literal slice, load the source word, shift
right and mask with 262143. When the slice
crosses a boundary, first load the second word,
shift left by the complementary offset and OR
before the mask. Crossing is determined at code
generation by (18*l mod 64)+18>64, giving three
crossings. The last slice fits within the padded
word. Multiply the key by k and add the literal
slice offset and direction base. Then for each
entry word add its literal word offset, load and
XOR into its accumulator. Store all accumulators
to the current cube array. The update costs

    Aupdate=2+2k+12*(6+3k)+3*3.

It is followed by flipping the bit of the field,
so the cube invariant holds for the child. Before
an inner backtrack restore the saved cube from the
row at current depth. Its cost is the same as the
save. Each edge therefore charges Aupdate plus
Asave=2*(2+3k).

Use the following template for the nine subset
trees of an orientation, with its own depth, cap,
length L and arrays. SLOT converts a logical slot
to its actual coordinate. It is a free move in
the anchor tree, and adds the list base and loads
one entry in an inner tree.

    VISIT:
        run the stated visit routine
    EXTEND:
        stop := eq(depth,cap)
        branch if stop to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        stop := eq(i,L)
        branch if stop to BACK
        j := add(i,1)
        store(adr,j)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        coordinate := SLOT(i)
        [non-anchor only] SAVE_AND_UPDATE(coordinate)
        FLIP(coordinate)
        jump to VISIT
    BACK:
        done := eq(depth,0)
        branch if done to DONE
        adr := add(LAST,depth)
        i := load(adr)
        [non-anchor only] RESTORE()
        depth := sub(depth,1)
        coordinate := SLOT(i)
        FLIP(coordinate)
        jump to EXTEND

FLIP shifts the actual coordinate right by 6 and
ANDs it with 63 for word and bit indices. Shift 1
left by the bit index, add the literal field array
base to the word index, load the word, XOR and
store. This costs seven, or at most nine with
SLOT. Fields keep zero padding since coordinates
are below 200.

A node advances its NEXT before its child, and the
child begins at the next larger slot. LAST records
the undo slot. Thus the tree visits exactly one
node per increasing subset up to its cap and undoes
each descent on return. For N visits the template's
event counts are N-1 successful extensions, N failed
extensions, N-1 nonroot returns and one root return.
They cost at most 23, 6, 15 and 2 respectively,
before the affine and visit macros. Their total is
at most 44N.

## Full primitive traversal, scoring and output

Extract Y_0..Y_4,P_0..P_4 into their padded arrays.
For each literal piece load and shift right from
its within-word input offset when nonzero. On a
crossing also load the next word, shift left by
the complementary offset and OR. Mask when the
piece width is below 64, then store to its padded
word. The 40 pieces use at most 280 primitives.
Load all received words from these arrays,
popcount and add to BestD, starting at zero, at
cost 120. This is the distance of the all-zero
fallback. Precompute I_i=Y_i^3 for all five
systematic positions with the cube macros.

Each orientation has a separately emitted block
with literal anchor and child field indices. To
start the anchor, copy its Y to X with 2k
operations. Move its depth to zero, store zero
into its NEXT[0], and fall through to its VISIT.
The anchor tree has L=200, cap=R. At each anchor
VISIT compute C from the current X with C0
operations. Add the anchor depth to DEPTH_BASE
and load its budget in two operations. Jump to
the first child setup with that remaining budget.

At a child setup with its literal field index j
and known predecessor cube, form the two lists.
Move the list lengths d and g to zero. At each
literal coordinate load the corresponding words of
Y_j, P_j and C_(j-1), XOR twice, then shift right
and AND with 1 for the disagreement indicator v.
Select d or g as the write offset and the
appropriate literal list base using v. Add them,
store the coordinate, and add v to d and 1 xor v
to g. This uses 14 primitives per coordinate,
including both increments and the XOR. The
coordinate lists occupy distinct 200-word arrays
for this child level.

For child level l=1,2,3,4 let Reserve be the literal
constant 5*(4-l). With Rem the current remaining
budget use the following nine primitives.

    Base := add(d,5)
    Needed := add(Base,Reserve)
    Skip := ltu(Rem,Needed)
    branch if Skip to the parent's EXTEND
    Raw := sub(Rem,Base)
    Available := sub(Raw,Reserve)
    r := div(Available,6)
    large := ltu(g,r)
    r := select(large,g,r)

The results after the branch are
only used on the non-skipped path. The total is at
most nine on either path. The subtraction is
nonnegative on that path. These steps choose the
greatest allowed h, capped at the agreement-list
length. Raw is kept for this active level.

If the setup is accepted, copy Y_j to X_j and I_j
to C_j, costing 4k. Move the disagreement depth
to zero and store zero at its NEXT[0]. Fall
through to its visit template with L=cap=d. At each
disagreement VISIT move the agreement depth to
zero, store zero at its own NEXT[0], and jump
to its VISIT. The agreement tree has L=g, cap=r.
All its subtree updates use its distinct saved
rows. At an agreement VISIT of levels one through three,
multiply the current agreement depth by six
and subtract from Raw to form the next Rem, then
jump to the next child setup. These are three
primitives. At level four
perform the score routine below, with continuation
at agreement EXTEND. The agreement template's
DONE is its disagreement EXTEND. The disagreement
DONE is the parent agreement EXTEND or anchor
EXTEND for level one. These returns use the
template's already counted branches. The inner
trees restore a field and its cube on each return,
and distinct levels' arrays do not change the
parent values.

At a level-four agreement VISIT, move Dist to
zero. For each literal field i and word position,
load X_i and the predecessor cube. XOR them to
form the candidate parity. Load its two received
systematic and parity words. For each of the two
words XOR the candidate and received, popcount,
and add to Dist. This uses 11 primitives per
field word. Compare Dist with BestD by ltu and
branch to agreement EXTEND if it is zero. On a
strict improvement move Dist to BestD, load each
candidate field and store its padded words to
B, in 2*5*k operations total. Fall through to
agreement EXTEND. This retains the first minimum.
The strict update occurs at most 2000 times over
all orientations because BestD starts at at most
2000 and decreases by an integer each time.

The anchor DONE is the next orientation's setup
or the final packing label for the last one.
Pack the five B arrays as increasing disjoint
output intervals. For each piece load its word
and store directly if its output target starts at
a word boundary. Otherwise shift left, load the
output word, OR and store. On a crossing also
shift the original piece right by the complementary
offset, load the next word, OR and store.
This uses at most nine primitives per piece, at
most 180 overall. All 1000 bits are written and
the top output padding stays zero. At aligned
starts no earlier piece has used that word. The
current cube accumulators use at most 2k field
registers. At each of four active child levels, save both
depths, the list lengths, cap and remaining budget
state in at most seven scalars. Besides those and
the 2k cube registers, the anchor depth, incumbent
score and all current instruction temporaries fit
in at most twenty additional registers. Field arrays,
lists, saved values and incumbent words are all
in memory, so the register bound is within 64.

## Worst-case count

For fixed received words and a fixed outer candidate
consider the raw tree of ALL 2^200 possible
assignments at each child, before budgets. Give
each assignment a conditional weight

    2^(-d) * (1/65)^h * (64/65)^(g-h),

where d is the current number of disagreements,
g=200-d and h the assignment's number of agreement
flips. For a fixed parent these weights sum to one
by binomial expansion over the g agreeing positions
and the 2^d disagreement assignments. They are a
counting measure on continuations for the fixed
received data. The product of such weights over
a prefix is its weight in a normalized measure at
that level, since every node's conditional weights
sum to one.

The checker verifies 65^200<=2^5*64^200.
Each child weight is at least 2^(-(5+d+6*h)):
the factor 64^(-h) equals 2^(-6*h) and
(64/65)^g is at least 2^(-5) since g<=200.
A visited prefix of length l satisfies total cost
at most K-5*(4-l) by the decoder's reserve.
It therefore has weight at least 2^(-(K-5*(4-l))).
The number of agreement visits at that level over
all setups per anchor is at most

    N_l = 2^(K-5*(4-l)).

Each such visit represents one distinct prefix
because the subset search enumerates each child
once with those parents. Put N_0=1. The number
of setups at level l is at most N_(l-1), because
each previous agreement visit can call it once.

A child setup excluding its subtree costs at most
E=14*200+9+4k+1. At each accepted setup there are
D disagreement visits and A agreement visits with
D<=A. Across these trees the number of successful
edges is (D-1)+(A-D)=A-1. A successful edge charges
one affine update and save/restore. The template
cost outside this is up to 46 per disagreement visit
including reset and jump, and 44 per agreement visit.
The additional agreement cost is three to pass the
remaining budget at the first three levels, or
11*5*k+2 for the final score outside the rare winning
updates. Charging D<=A permits the bounds

    Vmid = 44+46+3+Aupdate+Asave,
    Vleaf = 44+46+11*5*k+2+Aupdate+Asave.

Thus the cost per anchor with budget K is at most

    C0+3+44 + sum_(l=1..4)
        (E*N_(l-1)+V_l*N_l),

where V_1=V_2=V_3=Vmid and V_4=Vleaf. An orientation
begins with 2k+1 primitives, then sums the bound
over its anchor subsets grouped by their bands.
The five initial cubes cost 5*C0. The fixed cost
100000 covers extraction, fallback score, winning
updates and output. Summing all five orientations
gives

    W=980053103610888392142782257636072437922485562425655103461054793071585275413364993549286810492096915061475.

The verifier checks W^500<=2^174387 in integers.
All runtime trees are finite and the counted
schedule matches the stated deterministic decoder
on every possible received word.
