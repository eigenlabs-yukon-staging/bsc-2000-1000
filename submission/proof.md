# Cyclic cube parities with an integer inner budget

## Code and decoder

Use four 250-bit information fields b_0 through b_3
in little-endian order, with coordinates in the
polynomial basis of K=F2[X]/(P), where

    P = X^250 + X^59 + X^40 + X^32 + 1.

The declared computation verifies X^(2^250)=X modulo
P and gcd one between P and X^(2^k)+X for k=50,125.
An irreducible factor of degree d would have d divide
250 because its generating root is fixed by the 250th
Frobenius power. Any proper divisor of 250 divides
one of 50 and 125. Thus the two gcds prove irreducibility.

Transmit all four b_i followed by their four parity
fields in the same index order, with

    a_i = b_i xor b_(i-1 mod 4)^3.

This is an injection carrying exactly 1000 bits in
2000 binary channel uses. All b_i are independent and
uniform under uniform information.

The decoder begins with the all-zero candidate.
Write Y_i and P_i for the received systematic and
parity fields. Choose each i in turn as an anchor,
in order 0,1,2,3. Traverse all subsets of up to
R=31 flips to Y_i. For anchor weight u the
following inclusive bands determine the total
remaining integer budget K(u), in bits.

    0 <= u <= 16: K(u) = 295
    17 <= u <= 18: K(u) = 287
    19 <= u <= 20: K(u) = 282
    21 <= u <= 22: K(u) = 275
    23 <= u <= 24: K(u) = 269
    25 <= u <= 26: K(u) = 263
    27 <= u <= 28: K(u) = 257
    29 <= u <= 31: K(u) = 249

For each anchor candidate visit its other three
fields in cyclic order. At a child field j the
predecessor candidate x_(j-1) is known. Compare
Y_j with P_j xor x_(j-1)^3. Suppose they disagree
at d positions. Any proposed child x_j is specified
by flips relative to Y_j. If h flips are on agreeing
positions, the child costs

    6+d+6*h

units. Keep prefixes whose total cost is at most the
anchor budget minus six per not-yet-chosen field.
For each child this permits all disagreement flips
and all agreement subsets through a specific cap.
The reserve ensures possible remaining baseline
costs fit, so it does not exclude any complete path
within budget. A path of three children within the
total budget is scored by its full codeword distance.
Return the earliest minimum over both the fallback
and the four orientation searches. The precise order
is the subset tree order in the scalar program.

## Sufficient admitted shells

Conditional on a true predecessor, the two received
observations of a child differ from its true field
by the noises in its own systematic and parity fields.
If those noises have d disagreements and h double
flips, their group contributes weight d+2h. Given
systematic weight u, the corresponding number of
patterns is

    binom(250,d)*binom(250-d,h)*binom(d,u-h).

Along the true path a disagreement has no constraint,
and h is the number of needed agreement flips.
Thus its three children have exact total cost
18 plus the sum of their values d+6*h. The
reserve in the decoder permits those children if
this total is within the anchor's budget.

For the lower count choose a preferred true anchor
for each noise pattern as follows. Find the first
band occupied by any of the four systematic noise
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
specified band and index, and all four possible
preferred indices give disjoint events. The earlier
and later groups are the same three children of that
anchor, in some cyclic order. For actual true paths
their total cost is independent of the traversal order.

Make two polynomials of one child group in variables
(C,H). In each polynomial the monomial for d
disagreements and h double flips is C^(d+6*h)*H^h.
Its coefficient is the preceding binomial expression
summed over systematic weights at least low for
polynomial A, and greater than high for polynomial B.
For the four possible preferred-index placements,
the combined polynomial of the other three groups is

    A^3+A^2*B+A*B^2+B^3 = (A+B)*(A^2+B^2).

The budget truncates its C degree to K-18. For a
term C^c*H^h, d=c-6*h and its total noise weight
is d+2*h. Combine these child counts by ordinary
weight convolution with the preferred anchor's own
250 systematic and 250 unrestricted parity bits.
The resulting 231 coefficients count exactly the
sufficiently admitted patterns of total weight
0..230 for that preferred band. The eight bands'
sufficient events are disjoint.

In the experiment use stride s=1+floor((K-18)/6)
and encode C^c*H^h at degree s*c+h as a packed
univariate coefficient. Any term with h>=s has
c>=6h>K-18, so cannot alias a retained degree.
Terms with c>K-18 lie beyond all retained degrees.
For each term with c<=K-18, its true noise weight
is at most c. Therefore every retained coefficient
in a product of at most three group factors or in
the displayed sum is at most
4*binom(1500,K-18). This follows by
dominating by all unrestricted noise patterns of
that weight and using monotonicity of the binomial
coefficient in the indicated prefix. The declared
ranges have K-18<=750. The checker verifies
this bound is below the coefficient radix 2^1120.
Integer multiplication therefore computes these
truncated polynomials exactly.

The four declared runs each reproduce 231 coefficients
for two consecutive bands. Sector j checks bands
2*j and 2*j+1 against shells[j], the two arrays in
the fixed file shellsj.json. The final run sums all
eight independently checked arrays and computes
the full certificate. Define E_t to be these sufficiently
admitted noise patterns through weight t, and J_t
their cumulative number. The subsequent argument
only needs E_t to imply actual decoder admission.

## Competitors changing all four fields

Write q=2^250 and Q=q^4. For nonzero e, the cube
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
directions e_i in all four fields. Its parity
difference in field i is e_i plus the cube
difference in field i-1. Thus, for fixed directions,
the four translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 1000
parity bits can be evaluated by expanding the
four affine indicators. Multiplication by their
common point mass 16/Q leaves 1/Q times a sum
of sixteen terms. The empty mask term is the
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

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(250*h,w).

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
most binom(250*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(250*(4-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(250*(4-k),b) *
    H_k(t-a-b, binom(250*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(4,k) selections of masks, and add
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

A nonempty proper subset of changed systematic
fields on the cycle is a disjoint union of maximal
runs of consecutive positions, each of length a
in {1,2,3}. A run changes its a systematic fields
and a+1 parity fields: its own parity fields and
the next parity after the run. For different runs
these affected coordinates are disjoint. A
competitor differing on more than one run and no
further from y must have at least one run whose
local distance is no greater than the true word's.
The competitor which differs only on that run has
the same local distance and no change outside it,
so it is also no further from y. It suffices to
bound conflicts with a single run of changes.

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

    P_(a,d)(X) = 2^d*X^d*(1+X^2)^(250-d)*(1+X)^((2a-1)*250).

The repeated pair gives the first three factors.
The other (2a-1)*250 positions comprise the
remaining systematic and varying parity fields.
This polynomial counts all candidate patterns on
the pair, even ones implying a zero direction,
so is an upper envelope. With true local weight
w, the expected number of run competitors no
further than truth is at most 2^a/q^a times the
cumulative coefficients of this polynomial through w.
Clipping at one bounds the probability of a conflict
uniformly for each such noise. The unchanged fields
do not affect the distance comparison.

For each w construct S_(a,w) from this last bound
by multiplying its numerator min(q^a,2^a*volume)
by binom(250,d)*[X^w]P_(a,d) and summing over d.
There are eight code fields in all. Therefore

    T_(a,v) = sum_(w=0..v) S_(a,w)*binom(2000-(2a+1)*250,v-w)

counts an upper numerator for the full noise weight
v, including any patterns on untouched coordinates.
For a fixed run start and length a the probability
of a conflict at weight v is bounded by

    T_(a,v)*9^(2000-v)/(q^a*10^2000).

The verifier computes the unshifted coefficients
(1+X^2)^250*(1+X)^((2a-1)*250) and successively
divides by 1+X^2 as d increases, using the recurrence
new[i]=old[i]-new[i-2]. It uses the exact polynomial
prefixes through weight 230. There are
four possible run starts, so the union bound over all
lengths uses four times the sum of these bounds.

## Occupancy after removal

Among the admitted pairs (u,y) with noise in E_t,
remove any candidate for which there exists a
run competitor no further from y. Let L_t be the
cumulative number removed. If two distinct
retained candidates at a common y differed on a
proper subset of fields, take the one with greater
or equal distance. The other is a proper-change
competitor no further than it and some single run
would force removal, a contradiction. Thus all
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
      + 4 * sum_(a=1..3) sum_(t<=230) alpha_t*T_(a,t)*w_t/(q^a*10^2000),

where U_(-1)=0. Omitted true-success terms are
nonnegative. The declared computation checks this
quantity is at most 6975456245706939/125000000000000000 using exact arithmetic.

## RAM representation and field arithmetic

Each field is a padded little-endian array of k=4
words with h=58 live bits in the last word.
Construction.py specifies every address. Starting
at word 48 put the six groups of four field arrays
Y, P, X, C, I, B, in that order, then SQ. Y and P
store the received systematic and parity fields.
X and C hold the current candidates and their
cubes. I holds the cubes of the corresponding
received systematic fields. B holds the incumbent
message, initially zero. Each of seven trees has
a disjoint NEXT and LAST array of 251 words. After
those put the two 250-word coordinate lists for
each of the three child levels, then the two
251-by-k saved-cube arrays for each of those levels.
The multiplication scratch holds sixteen k-word
rows. Then come the fixed reduction, square, depth budget and affine tables. All are consecutive with
exclusive end 3670060068, within the 32-bit address
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
of c and P+X^250. This fits in one word since the
largest possible degree is 62. To multiply a field
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
c*X^250 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=32 byte slots, emit a
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
Move R registers to zero. For each of n0=63 nibble
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

For direction i the fixed affine slices cover 14
successive chunks of 18 bits, up to position 251,
with positions 250 and 251 set to zero. The entry
d of slice l XORs the linear columns
e*X^(2j)+e^2*X^j for all its selected live bits j.
It includes e^3 exactly for slice zero. Entry
word h0 is at

    AFFINE_BASE+i*AFFINE_STRIDE+l*262144*k+d*k+h0,

with AFFINE_STRIDE=14*262144*k. The generator
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

    Aupdate=2+2k+14*(6+3k)+3*3.

It is followed by flipping the bit of the field,
so the cube invariant holds for the child. Before
an inner backtrack restore the saved cube from the
row at current depth. Its cost is the same as the
save. Each edge therefore charges Aupdate plus
Asave=2*(2+3k).

Use the following template for the seven subset
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
are below 250.

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

Extract Y_0..Y_3,P_0..P_3 into their padded arrays.
For each literal piece load and shift right from
its within-word input offset when nonzero. On a
crossing also load the next word, shift left by
the complementary offset and OR. Mask when the
piece width is below 64, then store to its padded
word. The 32 pieces use at most 224 primitives.
Load all received words from these arrays,
popcount and add to BestD, starting at zero, at
cost 96. This is the distance of the all-zero
fallback. Precompute I_i=Y_i^3 for all four
systematic positions with the cube macros.

Each orientation has a separately emitted block
with literal anchor and child field indices. To
start the anchor, copy its Y to X with 2k
operations. Move its depth to zero, store zero
into its NEXT[0], and fall through to its VISIT.
The anchor tree has L=250, cap=R. At each anchor
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
coordinate lists occupy distinct 250-word arrays
for this child level.

For child level l=1,2,3 let Reserve be the literal
constant 6*(3-l). With Rem the current remaining
budget use the following nine primitives.

    Base := add(d,6)
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
rows. At an agreement VISIT of levels one or two,
multiply the current agreement depth by six
and subtract from Raw to form the next Rem, then
jump to the next child setup. These are three
primitives. At level three
perform the score routine below, with continuation
at agreement EXTEND. The agreement template's
DONE is its disagreement EXTEND. The disagreement
DONE is the parent agreement EXTEND or anchor
EXTEND for level one. These returns use the
template's already counted branches. The inner
trees restore a field and its cube on each return,
and distinct levels' arrays do not change the
parent values.

At a level-three agreement VISIT, move Dist to
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
B, in 2*4*k operations total. Fall through to
agreement EXTEND. This retains the first minimum.
The strict update occurs at most 2000 times over
all orientations because BestD starts at at most
2000 and decreases by an integer each time.

The anchor DONE is the next orientation's setup
or the final packing label for the last one.
Pack the four B arrays as increasing disjoint
output intervals. For each piece load its word
and store directly if its output target starts at
a word boundary. Otherwise shift left, load the
output word, OR and store. On a crossing also
shift the original piece right by the complementary
offset, load the next word, OR and store.
This uses at most nine primitives per piece, at
most 144 overall. All 1000 bits are written and
the top output padding stays zero. At aligned
starts no earlier piece has used that word. The
current cube accumulators use at most 2k field
registers. The seven tree depths and limits and
the remaining tests, counters and temporaries
need at most 44 further scalars. Field arrays,
lists, saved values and incumbent words are all
in memory, so the register bound is within 64.

## Worst-case count

For fixed received words and a fixed outer candidate
consider the raw tree of ALL 2^250 possible
assignments at each child, before budgets. Give
each assignment a conditional weight

    2^(-d) * (1/65)^h * (64/65)^(g-h),

where d is the current number of disagreements,
g=250-d and h the assignment's number of agreement
flips. For a fixed parent these weights sum to one
by binomial expansion over the g agreeing positions
and the 2^d disagreement assignments. They are a
counting measure on continuations for the fixed
received data. The product of such weights over
a prefix is its weight in a normalized measure at
that level, since every node's conditional weights
sum to one.

The checker verifies 65^250<=2^6*64^250.
Each child weight is at least 2^(-(6+d+6*h)):
the factor 64^(-h) equals 2^(-6*h) and
(64/65)^g is at least 2^(-6) since g<=250.
A visited prefix of length l satisfies total cost
at most K-6*(3-l) by the decoder's reserve.
It therefore has weight at least 2^(-(K-6*(3-l))).
The number of agreement visits at that level over
all setups per anchor is at most

    N_l = 2^(K-6*(3-l)).

Each such visit represents one distinct prefix
because the subset search enumerates each child
once with those parents. Put N_0=1. The number
of setups at level l is at most N_(l-1), because
each previous agreement visit can call it once.

A child setup excluding its subtree costs at most
E=14*250+9+4k+1. At each accepted setup there are
D disagreement visits and A agreement visits with
D<=A. Across these trees the number of successful
edges is (D-1)+(A-D)=A-1. A successful edge charges
one affine update and save/restore. The template
cost outside this is up to 46 per disagreement visit
including reset and jump, and 44 per agreement visit.
The additional agreement cost is three to pass the
remaining budget at the first two levels, or
11*4*k+2 for the final score outside the rare winning
updates. Charging D<=A permits the bounds

    Vmid = 44+46+3+Aupdate+Asave,
    Vleaf = 44+46+11*4*k+2+Aupdate+Asave.

Thus the cost per anchor with budget K is at most

    C0+3+44 + sum_(l=1..3)
        (E*N_(l-1)+V_l*N_l),

where V_1=V_2=Vmid and V_3=Vleaf. An orientation
begins with 2k+1 primitives, then sums the bound
over its anchor subsets grouped by their bands.
The four initial cubes cost 4*C0. The fixed cost
100000 covers extraction, fallback score, winning
updates and output. Summing all four orientations
gives

    W=35901097407128634856584690105857357386056636918659727367979507489736241899473650408828013846706605362461806433691848376.

The verifier checks W^125<=2^49229 in integers.
All runtime trees are finite and the counted
schedule matches the stated deterministic decoder
on every possible received word.
