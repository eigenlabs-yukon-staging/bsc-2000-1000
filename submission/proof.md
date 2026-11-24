# Coupled field cubes with two bounded searches

## Code and candidate set

Use little-endian bit order. Represent a 500-bit word as a
polynomial over F2 with coefficient of X^i at bit i. Work in
K=F2[X]/(P), where P=X^500+X^21+X^13+X^8+1. The declared
experiment verifies X^(2^500)=X modulo P and gcd one between
P and X^(2^k)+X for k=100 and 250. Each irreducible factor's
degree must divide 500, since its root generates its field
and is fixed by the 500th Frobenius power. Every proper degree
dividing 500 divides 100 or 250 because the prime divisors of
500 are 2 and 5. The gcd tests therefore prove irreducibility.

Let b,c be the low and high 500-bit halves of the information
word. Transmit four consecutive 500-bit fields indexed as
(a,b,c,d), with

    a = c xor b^3,     d = b xor c^3.

This carries the 1000 bits in exactly 2000 uses and is injective
because its middle two fields expose the original information.

Write (ya,yb,yc,yd) for the received fields. The decoder starts
with the valid candidate b=c=0. Its first search enumerates
every b at Hamming distance at most R=64 from yb. For an
outer candidate b, measure the distance i from yb by the
number of flips in its enumerated subset. Look up the budget
K(i) given by the following inclusive bands.

    0 <= i <= 52: K(i) = 181
    53 <= i <= 56: K(i) = 169
    57 <= i <= 60: K(i) = 157
    61 <= i <= 64: K(i) = 145

Compare yc with ya xor b^3. Their j disagreements divide the
500 inner coordinates into a set D of size j and a set G of
size 500-j. If j>K(i), this outer visit adds no candidates.
Otherwise choose the largest r in [0,500-j] such that

    2^j * sum_(h=0..r) binom(500-j,h) <= 2^K(i).

Starting the inner candidate at yc, enumerate all choices of
flips on D and every subset of at most r flips on G. Score
each resulting (b,c) by the full Hamming distance of its
codeword to y. The second search interchanges b and c and
interchanges a and d. It uses exactly the same bands and
inner rule, comparing yb to yd xor c^3 for each outer c.

Return the earliest candidate at smallest full distance,
with the zero candidate first, the b-outer search next, and
the c-outer search last. The primitive program below specifies
the deterministic order within each search. Both searches
are finite and the fallback exists on every received word.

## Counting admitted noise

For an actual message, if its b field is among the first
search's outer candidates, the two observations for c differ
from c by the independent noises of the c and a fields.
For their two noise masks let j be their number of
disagreements and h their number of pairs with two flips.
The true inner word is among this outer visit's candidates
exactly when j is not skipped and h<=r(j,K(i)). There is no
condition on the first or second flipped bit at a disagreement
because every assignment there is covered. The same argument
applies to the b,d noise pair in the other orientation.

Use two 1000-bit noise groups, each consisting of a systematic
field and its partner in this rule. These groups are (b,d) and
(c,a). Let w be a group's total weight and u its systematic
weight. Given j disagreements and h double flips, w=j+2h.
The number of patterns with prescribed u is

    binom(500,j) * binom(500-j,h) * binom(j,u-h).

The final factor chooses which disagreements put their flip
in the systematic word. A band indexed by e specifies a
range of u and its budget K_e. Let O_e(w) count all group
patterns of weight w with u in that band. Let I_e(w) count
all weight-w patterns passing the inner cap at budget K_e,
with no condition on u. Let H_(e,f)(w) count group patterns
with u in band e that pass the inner cap at budget K_f.
These are one-dimensional integer arrays, evaluated from the
formula above or ordinary binomial convolution for O.

The count of weight-t patterns admitting the true message
in the b-outer search is the coefficient of degree t in
sum_e O_e*I_e, with multiplication denoting ordinary
convolution. It is the same for the c-outer search by
symmetry. Their intersection has coefficient from
sum_(e,f) H_(e,f)*H_(f,e). The two outer weights determine
their bands uniquely, which justifies summing the disjoint
intersection cases. Inclusion and exclusion gives the
exact admitted weight counts A(t). Define E_t to consist of
admitted noises of weight at most t, and J_t=sum_(v<=t) A(v).
All of these definitions depend only on the noise pattern.
For any u and z in E_t the decoder visits u on receipt
f(u) xor z, so the set of actual candidates contains every
pair counted here.

## Double changes and affine parity distributions

Write q=2^500 for the field size and Q=q^2 for the number
of messages. The pair (b,c) is independent uniform in K^2.
For a nonzero field difference e, the cube difference is

    (x+e)^3 + x^3 = e^3*((x/e)^2 + x/e + 1).

The map v->v^2+v is F2-linear with kernel of size two.
Its image is the kernel of absolute trace. The trace is
nonzero because its polynomial has degree 2^499<q and
leading coefficient one, and its values are in F2. The
trace vanishes on v^2+v by the Frobenius identity. Thus
the cube difference is uniform on an affine trace
hyperplane, with each point having mass 2/q.

The hyperplane's normal in the binary word basis is the
nonzero mask with entries Tr(e^-3 X^j). The map taking a
field multiplier to that coordinate vector is injective.
A nonzero multiplier permutes K, so could give zero on
every basis vector only if trace vanished on all K.
Consequently any mask has at most three original nonzero
directions, since e^-3 fixed is a cubic equation in e.

For a competitor (b+p,c+r) with p,r both nonzero, its
parity differences at a and d are respectively r plus
the cube difference at b with direction p, and p plus
the cube difference at c with direction r. Conditional
on the fixed directions these two affine distributions
are independent, because b and c are independent. Adding
their known true parity noises merely translates the two
hyperplanes. For a 1000-bit parity ball of radius l, the
expected indicator of membership in the ball is at most

    (V_l + abs(K_l(w_p)) + abs(K_l(w_r))
       + abs(K_l(w_p+w_r))) / Q,

where V_l is its volume and w_p,w_r are the weights of
the 500-bit normal masks. To see this, expand the product
of the two affine hyperplane indicators, using
(1+sign*(-1)^(mask dot x))/2 for each one. Multiply its
intersection size by 4/Q. The result is the ball volume
and three signed terms, all bounded by the displayed
absolute magnitudes. A 1000-bit mask of weight w has
signed ball count

    K_l(w) = [X^l] (1-X)^(w-1)*(1+X)^(1000-w).

This follows by summing signed shell coefficients through
l, which divides the shell polynomial by 1-X. The exact
coefficient recurrence follows by differentiating this
polynomial F. With N=999 and g=N-2*(w-1) it satisfies
(1-X^2)F'=(g-NX)F, hence

    (l+1)*K_(l+1) = g*K_l - (N-l+1)*K_(l-1),

with K_(-1)=0 and K_0=1. The computation performs exact
integer division and checks every division it uses.

For a fixed l let H1(l,n) be the sum of the largest n
magnitudes in the one-half mask list, with multiplicity
3*binom(500,w) at each w>0. The sum of the largest n
two-half-mask magnitudes is H2(l,n), where the weight-w
multiplicity is

    9 * (binom(1000,w) - 2*binom(500,w)),   w>0.

This counts pairs of nonzero masks of combined weight w
with up to nine ordered direction pairs per pair of masks.
The subtracted terms are the cases where one half is zero.
Both lists are evaluated by a descending sort of groups,
using cumulative counts and magnitudes and a partial last
group. Their multiplicities exceed the capacities used.

Now fix a true noise z of total weight at most t. A
double-change competitor's systematic errors have weights
i=wt(z_b+p) and j=wt(z_c+r). There are at most n_i and n_j
possible directions on the two sides, with n_i=binom(500,i).
For i+j<=t its available parity-ball radius is l=t-i-j.
The sum of the three imbalance terms over directions
in this slice is bounded by

    n_i*H1(l,n_j) + n_j*H1(l,n_i) + H2(l,n_i*n_j).

The zeroth term contributes n_i*n_j*V_l. For a
competitor arising from E_t, either its b or its c
systematic weight is at most R, since its error is
admitted in at least one orientation. Thus it suffices
to sum the zeroth and imbalance terms over
i+j<=t and min(i,j)<=R. Call the resulting integer
D_t. This bounds the expected number of double-change
competitors at distance at most t and with admitted
error by D_t/Q, for fixed z and uniform true
message. The calculation checks that D_t is
nondecreasing through the cutoff t=230.

## Single changes and probability of removal

For a change in c alone, the a and c fields change by
the same nonzero direction r, and b is fixed. The
remaining parity field d changes according to the cube
difference at uniform c. Thus its error is an affine
hyperplane distribution of mass 2/q per point. The
single direction in b alone has the same bounds by
symmetry.

Fix the true noise in the a,c,d fields, and let j be
the number of disagreeing positions in the paired a,c
fields. Applying a common difference to them preserves
their set of disagreements. On this fixed set there are
two patterns of weight one per position. Outside it
there is one pair of weight zero and one of weight two.
Including the remaining d field, the weight polynomial
for patterns with the given disagreement positions is

    P_j(X) = 2^j * X^j * (1+X^2)^(500-j) * (1+X)^500.

A choice of changed errors on the two paired fields
determines the common direction uniquely. For a full
three-field true weight w, the expected number of
single-change competitors no further than the truth is
bounded by 2/q times the cumulative coefficients of P_j
through w. This uses the point-mass bound on d. The
unchanged b field adds the same cost to both words.
By the union bound and clipping at one, the conditional
failure over the uniform true message is bounded by

    min(q,2*sum_(h<=w) [X^h]P_j)/q

for each fixed three-field noise pattern. For a fixed
three-field weight w, multiply the numerator by
binom(500,j)*[X^w]P_j and sum over j to get S_w.
In particular the sum, over all three-field errors
of weight w, of their respective conditional
failure bounds is at most S_w/q.

With the additional unchanged b noise field, define

    T_v = sum_(w=0..v) S_w*binom(500,v-w).

The upper probability of this type of single-change
conflict and total true noise weight v is
T_v*9^(2000-v)/(q*10^2000). The convolution simply
accounts for all patterns on the unchanged systematic
field, which is independent of the three other noises
and whose weight does not change the pair comparison.
Both single-change types together have probability at
most twice this bound for any given v.

The computation forms the base coefficients of
(1+X^2)^500*(1+X)^500. When j increases one, divide
the unshifted polynomial by 1+X^2 on the truncated
coefficients using new[k]=old[k]-new[k-2], with
negative indices zero. This generates the exact
terms of P_j and only weights through v=230 are
needed.

## Occupancy and weight-dependent removal

Fix a total radius t. For each received y, consider the
candidate pairs (u,y) with noise in E_t. The total number
over y is Q*J_t. For this argument remove every such
pair whose true message has any competing codeword
sharing b or c that is no further from y. Let L_t
count the removed pairs cumulatively through radius t.
Between two distinct retained messages at the same y
both coordinates must change. If they shared one,
the greater of their distances would fail retention,
with ties both removed. The retained competing
candidate itself has its error in E_t. Summing the
bounded double count D_t/Q over each of the J_t
noises and each of the Q true messages bounds the
total ordered distinct pairs of retained candidates
at any common y by J_t*D_t.

Let F_t be the number of received words with an
actual decoder candidate of distance at most t.
All received words with retained candidates are
included. Write m_y for the retained counts. For any
integer k>=1 and m>=0,

    k*(k+1)*1[m>0] >= 2*k*m-m*(m-1).

For positive m rewrite the inequality as
(m-k)*(m-k-1)>=0. At m=0 both sides vanish.
Use sum_y m_y=Q*J_t-L_t and the ordered pair bound
and divide the resulting inequality by Q. With
k_t=1+floor(D_t/Q), put

    G_t = floor(J_t*(2*k_t*Q-D_t)/(Q*k_t*(k_t+1)))
    alpha_t = 2/(k_t+1).

Then F_t/Q >= G_t-alpha_t*L_t/Q. The verified
monotonicity of D_t implies alpha_t is
nonincreasing. Let U_t=max_(i<=t) G_i and choose a
specific i(t)<=t attaining the maximum. Since F_t
is nondecreasing, it satisfies
F_t/Q >= U_t-alpha_(i(t))*L_(i(t))/Q.

On each received word the joint success probability
is the likelihood of its one winning candidate
divided by Q. The weight-t factor is 9^(2000-t)
divided by 10^2000. For the unnormalized factors
w_t=9^(2000-t) through t=230 set w_231=0. Summation
by parts weights each cumulative F_t/Q by the
nonnegative coefficient w_t-w_(t+1). Substituting
the bound from U and subtracting the loss gives
the following conservative treatment of a single
removed pair of true weight v. It occurs in
L_(i(t)) only for terms with i(t)>=v, hence t>=v
and alpha_(i(t))<=alpha_v. Its total contribution
to the loss numerator is at most alpha_v*w_v/Q
by telescoping the nonnegative coefficients over
t>=v. This argument applies also if no t attains
i(t)>=v. Dividing by the likelihood denominator
10^2000 and summing over removed pairs, the full
loss is at most

    2*sum_(v=0..230) alpha_v*T_v*w_v / (q*10^2000).

Here the two single-change directions use the bound
from the preceding section for each exact noise
weight. The nonnegative actual-success terms beyond
t=230 may be omitted. Consequently the exact
certified success lower bound is the sum of the
increments of U_t against w_t/10^2000 minus the
displayed loss. The experiment verifies that one
minus this lower bound is at most 5559788863928901/100000000000000000.

## Padded field words and fixed memory

Each field uses k=8 consecutive words in little-endian
order, with only h=52 live bits in the final word.
The layout defined by construction.py uses addresses
48 onward after the protected input and output prefix.
It has the twelve named k-word arrays in the stated
order, the three disjoint NEXT/LAST pairs each of
length 501 per array, two 500-word coordinate lists, two (501)-by-k
saved-cube arrays for the inner trees, sixteen
k-word multiplication rows, and the fixed reduction,
square, cap and affine difference tables. They
are consecutive, with exclusive memory end 2359439405.
All writable scratch starts at zero, has no aliases
with preloaded tables, and is overwritten on use.
The best b,c arrays start at all zeros for the fallback.

The two reduction tables for shifts by s=1 and s=4
have 2^s words each, with entry c equal to the carryless
product of c and P+X^500. The square tables, one per
byte position, have 256 rows of k words each. Byte
position j and byte d yield the field words of
(d*X^(8j))^2. The host data generator interleaves input
bits with zeros for a square and folds high polynomial
parts using the sparse modulus. A fold reduces degree
and preserves remainder. These entry definitions fix
every word of the finite table.

The depth table indexed by 0..R gives the base of
the cap row for that depth's band. Each row stores
r(j,K) for j=0..500, with a sentinel 501 for skipped
positions. Construction.py determines those values
using the exact cardinality inequality and the
experiment checks the corresponding integer budget
inequalities. No data-dependent table construction
occurs at decoding time.

## Primitive cube macros

All operations below use the bsc-ops-v1 primitives.
Literal loops and literal choices expand at code
generation. Moves of constants and registers are free.
A macro is expanded inline. There are no runtime
call or return instructions other than the explicit
branches. For a literal shift s in {1,4},
SHIFT_s acts on k field registers R0..R[k-1] as follows.

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

The descending order preserves the old lower words for
carries. The removed bits above coordinate 499 are
replaced by their modular contribution to the low
word. This contribution fits 64 bits because
deg(P+X^500)=21 and s<=4. The count is 3k+3.

To square a source array X0, initialize k register
accumulators to zero by moves. For each of the b0=63
literal byte positions, load the source word, shift
right to its byte, AND with 255, multiply the byte by
k and add the literal square table base for that
position. Then for each of its k field words add the
literal word index, load and XOR with the corresponding
accumulator. Store the accumulated value in SQ.
This costs S=b0*(5+3k)+k=1835. The map is correct
because squaring is F2-linear and the source high
padding supplies zeros to its final partial byte.

To multiply SQ by X0 into the specified cube array,
load SQ into k field registers and store the zero
vector and SQ into multiplication rows zero and one.
Emit SHIFT_1 three times to form rows 2,4,8, storing
each k-word row. For each of the other eleven nonzero
rows, in increasing order, load and XOR the rows for
its least set bit and remaining digits and store
the result. This costs four primitives per field word
for each of those eleven rows. Then move all R
registers to zero. For each of the n0=125 literal
nibble positions of X0 in decreasing order, load its
word, shift right to the nibble, AND with 15,
multiply by k and add the MULT base. Emit SHIFT_4
and XOR the indicated row into R, using one add,
load and XOR per field word. This is the Horner
evaluation of the product in the field. Store the
result to the specified cube array. The count is

    Mmul=3k+3*(3k+3+k)+44k+n0*(5+3k+3+3k)+k=7489.

The two macros together compute the exact cube. They
use no address outside the declared banks and always
leave padding zero.

## Incremental inner cubes

The two inner subset trees each maintain ICUBE as
the cube of INNER. Before a descent along actual
coordinate i, save all k words of ICUBE at row
depth of that tree's saved-cube array. Compute the
row base by mul and add and then use an add, load
and store per word, for 2+3k operations. Let e=X^i
and let x be the parent INNER. The cube change on
this edge is

    (x+e)^3+x^3=e*x^2+e^2*x+e^3.

The linear column at bit j is e*X^(2j)+e^2*X^j.
The affine table for i uses 36 slices of 14 input
bits, with its last four columns zero beyond the
field. The entry at slice l and value d XORs the
columns selected by d, and additionally includes
e^3 if l=0. It is stored as k words. The address of
its word h0 is

    AFFINE_BASE+i*AFFINE_STRIDE+l*16384*k+d*k+h0

where AFFINE_STRIDE=36*16384*k. The fixed entry
generator supplies exactly these values. Its
region follows the cap tables.

For the update compute the direction base with
mul and add, and load the k current ICUBE words
as accumulators. For each of the 36 literal slices
load its source word from INNER, shift right to
the slice start and AND with 16383. For a crossing
slice, also load the next word, shift it left by
the complementary offset and OR it into the value
before masking. The crossing indices satisfy
(14*l mod 64)+14>64, and there are six of them.
The final slice takes its extra zeros from the
padded field word. Multiply the key by k and add
the literal slice offset and the direction base.
For each output word add its literal displacement
to the entry address, load and XOR into the
corresponding accumulator. Store the k totals
to ICUBE. The update count is

    Aupdate=2+2k+36*(6+3k)+3*6=1116.

Then FLIP changes INNER at i, making ICUBE correct
for the child. On the backtrack RESTORE reads
the old saved-cube row, using mul and add and
an add, load and store per word. Restore before
lowering depth and flipping the same actual
coordinate to recover the parent. SAVE and
RESTORE together cost Asave=2*(2+3k)=52 for one
completed edge. No affine work is done by the
outer traversal.

## Complete search and score program

Initially extract each padded piece of ya,yb,yc,yd
from its literal input offset. Load its word and
shift right by the within-word offset when needed.
When the piece crosses a word boundary, also load
the next word, shift left by the complementary
offset and OR into the piece. AND with its width
mask if smaller than 64 and store to its Y array.
This costs at most seven primitives per piece, or
224. Initialize BestD to the distance of the zero
fallback. Load each Y word, popcount and add to
BestD, at cost 96. Keep BestD across both searches.

Emit separate b-outer and c-outer blocks with the
same schedule, switching the roles of the input
arrays as specified at the top. Start each block by
copying its systematic received outer field to
OUTER with 2k loads and stores total, moving its
outer depth to zero, and storing zero into its
outer NEXT[0]. Compute BASECUBE from the received
systematic inner field with the cube macros once
per orientation. Use the following subset tree
template wherever a bounded subset is needed.

    VISIT:
        emit the visit routine for this tree
    EXTEND:
        stop := eq(depth,cap)
        branch if stop to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        stop := eq(i,L)
        branch if stop to BACK
        next := add(i,1)
        store(adr,next)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,next)
        adr := add(LAST,depth)
        store(adr,i)
        coordinate := LOGICAL_TO_COORDINATE(i)
        [inner trees only] SAVE_AND_UPDATE(coordinate)
        FLIP(coordinate)
        jump to VISIT
    BACK:
        done := eq(depth,0)
        branch if done to DONE
        adr := add(LAST,depth)
        i := load(adr)
        [inner trees only] RESTORE()
        depth := sub(depth,1)
        coordinate := LOGICAL_TO_COORDINATE(i)
        FLIP(coordinate)
        jump to EXTEND

The length L is the number of available coordinate
slots, in the order defined for each tree. In the
outer tree take L=500, cap=R and the slots are
the coordinates themselves, and mapping a slot is a
free move. In the two inner trees the slot mapping
uses two instructions, adding the literal list base
and loading its i-th entry. FLIP shifts the actual
coordinate right by 6 for a word index, ANDs it with
63, shifts 1 left by that offset, adds the appropriate
field base to the word index, loads the word, XORs
and stores. The mapping and FLIP cost at most nine.
All coordinates are below 500. Each descent visits
one increasing subset and its return restores the
parent. Thus the template visits each subset of at
most cap positions exactly once, including the root.
For N visits count N-1 successful extensions of at
most 23 primitives, N failed attempts at most six,
N-1 nonroot returns of at most 15, and a root return
of two. Excluding the affine macros, the total is at most 44N.

At each outer VISIT, compute OCUBE from OUTER using
the cube macros. Take its two inner observations from
the received systematic inner field and the partner
field xor OCUBE. Initialize j and g to zero by moves.
For each literal coordinate i, emit three loads of
the relevant words of the systematic, partner and
OCUBE arrays, two XORs, a right shift and an AND of
one bit to obtain its disagreement flag v. Select j
or g by v for the write offset, select the literal
DISAGREE or AGREE base by v, add the chosen offset,
and store coordinate i. Increment j by v and g by
(1 xor v), charging the two adds and the XOR.
This uses 14 primitives per coordinate and gives both
disjoint sorted lists. Load the cap table row by
adding outer depth to DEPTH_BASE and loading, then
add j and load r. Compare with 501 and branch to
outer EXTEND on a skip. This totals six primitives.

Otherwise copy the systematic received inner field
to INNER and BASECUBE to ICUBE in 4k operations. Move the disagreement
depth to zero and store zero into its NEXT[0]. Fall
through to a disagreement subset tree with L=cap=j,
using the DISAGREE list for FLIP on INNER. Each of its
2^j VISITs moves the agreement depth to zero and stores
zero into its NEXT[0], then falls through to the
agreement subset tree with L=g and cap=r. Use the
AGREE list for its FLIP on INNER. Its DONE label is
the disagreement EXTEND, and the disagreement DONE
is the outer EXTEND. These targets use the existing
root-return branches. The inner array and its
cube are restored when leaving each subtree, since
both trees use the save and restoration macros with
separate saved rows. The exact number of
agreement VISITs at this outer candidate is
2^j*sum_(h<=r) binom(g,h), at most 2^K(i).
The number of disagreement VISITs is no larger.

At each agreement VISIT the maintained ICUBE is
exactly the cube of INNER. For the literal
chosen orientation OUTER and INNER represent b,c
in the appropriate order, and OCUBE and ICUBE
represent their corresponding cubes. Move Dist to
zero. For each literal field word position load
the four current field and cube words and XOR to
form the two parity words. Load the four received
words. For each coordinate array, take the XOR
difference, popcount and add it to Dist. This uses
22k primitives over all words. Compare Dist with
BestD by ltu and branch past the update if the result
is zero. Otherwise move Dist into BestD, load b,c
word by word and store them to BEST_B and BEST_C at
cost 4k. The score and strict update implement the
stated tie order. Because BestD is an integer
initialized at no more than 2000 and decreases only
on improvement, the 4k part runs at most 2000 times
over the entire decoder. All distances are exact,
including the zero padding.

The b-outer tree's DONE is the c-outer initialization.
The c-outer tree's DONE points to final packing.
Copy the two BEST fields to the output prefix in
order as disjoint intervals. For each piece load the
word and, if its target is aligned, store directly.
Otherwise shift left, load the previous output word,
OR and store. If it crosses, also shift the piece
right by the complementary offset, load the next
word, OR and store. There are at most nine primitives
per piece, or 144. In the aligned case no earlier
piece occupied that word. Initially all output words
are zero. The increasing disjoint bit ranges include
all 1000 bits and preserve zero padding in the output
word 47. The padded high words in BEST also remain
zero above their respective live bits. The field data arrays are in memory. At most 2k field
registers and 32 scalars are live at once, including
all three depths, limits, BestD, and temporaries.
This is within 64 registers.

Let N=sum_(i<=R) binom(500,i) and
V=sum_e 2^K_e*sum_(i in band e) binom(500,i).
Each outer search has N visits. Per outer visit
charge S+Mmul for its cube, 14*500 for both lists,
six for the table selection, 4k+1 for initialization
of the inner tree, and 44 for its traversal share.
Charge at most V agreement visits and V disagreement
visits in one orientation. Agreement work outside
affine operations is at most 22k+2+44 per visit.
Disagreement work outside affine operations costs
up to 45 including the reset store.

For a non-skipped outer node let D be its exact
number of disagreement visits and A its total
agreement visits, with A>=D. Its disagreement tree
has D-1 successful edges. Its D separate agreement
trees have A-D successful edges in total because
each tree has one root. Thus the sum of successful
inner edges is A-1, and each such edge has one
affine update and one save/restore pair. Over one
orientation charge at most V*(Aupdate+Asave).
The fixed cost 100000 covers all strict global
winning updates, input extraction, fallback score
and output. Each orientation adds initialization
2k+1 plus S+Mmul for BASECUBE. Thus

    W = 100000 + 2*(2k+1+S+Mmul)
        + 2*( (S+Mmul+14*500+6+4k+1+44)*N
             +(22k+2+44+45+Aupdate+Asave)*V )
      = 6231188497493191321731737558812847005717738060746248689479391340171171621792430053063998524982297192448920245941527907735028916776.

The exact exponent comparison is W^1000<=2^431169.
This is a worst-case count for the primitive decoder
on every received word, including words where a cap
causes one or more visits to skip.
