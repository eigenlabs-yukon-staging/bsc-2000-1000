# Three cyclic cube fields with two shortened coordinates

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

Each L consists of twenty-seven width-12 cells, two width-4 cells,
and two width-1 identity cells in increasing bit order. In a
width-12 cell with systematic value u, the low eleven parity bits
are the coefficients of the polynomial remainder of u*X^11 modulo
1+X^2+X^4+X^5+X^6+X^10+X^11. Its top parity bit is the sum
modulo two of the Hamming weights of u and the low eleven bits.
The width-4 cell uses columns 14,13,11,7 as binary vectors.
Consequently L is linear, and the high two coordinates of L(b0)
are zero.

The decoder visits each anchor x0 in U within 68 systematic
errors of Y0. At depth d (the anchor systematic distance), use
C(d)=4340-40d. For a fixed anchor, visit child fields 1 and 2
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
All possible bits are retained in each row and sorted by distance
and then by numerical error, so cost is nondecreasing.

There are 62 total child cells. Continue with a prefix of length
i only if its used cost is at most C(d)-(62-i). Every complete
tuple of cost at most C(d) passes all prefix tests, since every
later cell costs at least one. Compare all returned messages by
their full distance. Initialize the incumbent to zero, replace
only on a strict improvement, and use the traversal order below.
This defines the selected message on every received word.

## Exact admission and comparison with full-code ML

On the true child path the actual predecessor cubes are known
at each cell, so adjustment gives the original independent
systematic and parity noise. For a paired local syndrome define
h_s(w)=#{e:D_s(e)=w}. Then the numerator generating polynomial of
its true integer cost has a coefficient contribution
h_s(w)*9^(2m-w) at c_s(w), summed over all s,w. For an
unpaired cell the numerator has terms 9 at c(10,9), and 1
at c(10,1). Thus the child cost law has numerator equal to a
product of 54 width-12, 4 width-4, 2 paired width-1 and
2 unpaired factors, with denominator 10^1334. The experiment
computes every histogram by full enumeration and groups equal
rows with their exact multiplicities. Its packed digits use 4500
bits, exceeding the bit length of 10^1334. Each partial product
has a sum of coefficients at most this numerator bound, so
truncated integer multiplication gives exact coefficients.

The child noises are disjoint from the 332 anchor systematic
noise bits. For d at most 68, multiply the child cost prefix
through C(d) by binom(332,d)*9^(332-d)/10^332.
The sum is the true-path admission probability of the bounded
traversal. Its complement E_miss also bounds absence from the
candidate set including the zero fallback.

For each received word choose a full-code ML minimizer with
preference for the restricted decoder's output if it is a
full-code minimizer. The ML likelihood is independent of which
minimizer is chosen because messages are uniform. If this
comparator chooses the true message and the true message is in
the restricted candidate list, the restricted output is itself a
full-code minimizer and is therefore the same chosen true message.
The actual error is at most E_miss plus the full-code ML error.

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
therefore have thirty-one small factors in the first-pair polynomial H.

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

Form the product of the thirty-one polynomials so obtained and
(1+Z)^u. At each local weight w multiply by the rational
common factors above, dividing by the thirty-one powers of 2^64
for the roots. This bounds the sum of the individual probability
bounds over true local noise patterns of that weight. Minimize
the result at each w over the stated grid and the trivial
binom(n,w) bound. Each candidate is uniform for every fixed
outside noise of that weight class, so the minimum before outside
expansion is valid. The local sum is then multiplied by
binom(2000-n,v-w) and channel factor 9^(2000-v)/10^2000.
Sum over w to get a bound J_(i,a,v) on each joint probability
for start i, length a and total true weight v. The checker
groups exactly identical geometries by their multiplicities.

The root numerator is at most 2^71 since each cell has no
more than 2^12 terms at r^2<=1. Each product of paired
components sums to at most 2^(668+31*71); the free factor
adds at most 2^1002. Thus 4000-bit digits suffice throughout
the packed polynomial products. All further arithmetic is rational.

## Full-code occupancy

At distance cutoff (t\le230), remove the pairs of true word and
noise whose noise weight is at most (t) and for which a single-run competitor
has distance no greater. Write (L_t) for the number of
these pairs over all (Q) messages. Two retained candidates at a
common received word must then differ in all three fields: otherwise
the farther one would be removable, as would both on a tie. There are (Q V_t-L_t) remaining
pairs, where (V_t=\sum_{s\le t}\binom{2000}{s}).
By the full-change estimate, their number of ordered distinct
pairs at common received words is at most (V_tD_t), obtained by
summing (D_t/Q) over all (Q V_t) true pairs.

For each (y) let (m_y) count its retained pairs, and let (F_t)
count received words with a full-code candidate within distance
(t). For every integer (k\ge1),

    k(k+1) 1[m_y>0] >= 2k m_y - m_y(m_y-1),

because on positive integers the difference is
((m_y-k)(m_y-k-1)\ge0). Summing and using the counts above gives

    F_t/Q >= G_t - alpha_t L_t/Q
    k_t = 1 + floor(D_t/Q),  alpha_t = 2/(k_t+1)
    G_t = floor(V_t*(2*k_t*Q-D_t)/(Q*k_t*(k_t+1))).

Set (U_t=\max_{i\le t}G_i) and choose a maximizing index (i(t)\le t).
By monotonicity of (F_t), it suffices to subtract the charge
(\alpha_{i(t)}L_{i(t)}/Q) from (U_t). Set
(w_t=9^{2000-t}) up to (t=230) and (w_{231}=0).
The exact success of full ML through the cutoff equals
\(\sum_t (w_t-w_{t+1})F_t/(Q10^{2000})\), because each received word
is decoded to one minimizer and contributes that likelihood divided
by (Q). In this sum a removed pair of true weight (v) can
be charged only if (i(t)\ge v), so (t\ge v) and
\(\alpha_{i(t)}\le\alpha_v\) by the nondecreasing (D_t).
Its total charge is at most (w_v\alpha_v/(Q10^{2000})).
Counting each single run by the stated probability bound gives

    E_ML <= 1 - sum_t (U_t-U_(t-1))*w_t/10^2000
            + sum_(i=0..2) sum_(a=1..2) sum_(v<=230) alpha_v*J_(i,a,v),

with (U_{-1}=0). The omitted contribution to success beyond
the cutoff is nonnegative. Both bounds used in the final error
score are checked with integer and rational arithmetic.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The six received fields are copied to padded little-endian
six-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, five frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed tables then contain reduction, square, local parity,
and sorted candidate rows. The Golay table has 2^24 words,
the Hamming table has 2^8 words, the repetition table has
2^2 words, and the unpaired systematic row has two words. Each word is an error in bits 0..11,
a distance in bits 12..16, and its integer cost in higher bits.
The generator specifies all table addresses within the 32-bit space.

Ranges with literal indices below expand when emitting the program.
References to literal RAM cells cost one load or store. Registers
and constants can be moved or substituted freely. All computed
addresses use the indicated chargeable arithmetic. Each conditional
transfer includes a comparison when one is shown and a branch.
At the boundary of each expanded block, other than the four
persistent registers (anchor depth, incumbent distance, remaining cost
and current distance), intermediate values
can be discarded. All 62 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use six accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=6 words and h=14 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^334. This fits in one word since the
largest possible degree is 8. To multiply a field
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
c*X^334 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=42 byte slots, emit a
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
Move R registers to zero. For each of n0=84 nibble
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
Use this for the 36 padded words of Y and P, using the exact
transmitted widths to place the field boundaries, then store each.
The shortened Y0 and P1 use two additional zero padding bits.
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
at most 3*(6+2)+1 primitives per word. Use the message field lengths 332,334,334. The final partial
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
local parity table base, load the parity image, and pack
these images into a six-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining=C(Depth) with one multiplication by 40 and one
subtraction from 4340. Move Distance=Depth and jump to CELL_INIT.
This uses one cube and at most 410 further primitives:
the 31 cells use at most 31*(6+2+4), and finishing AP
uses 18, two budget operations, and one jump.

### Cell traversal

Emit one literal block for each of the 62 cells, in
field order, then increasing offset within each field.
For this block let DepthCell be its 1-based cell index
and let m, j and i be its width, field and within-field
bit offset. Use frame RAM slots Base, U, Rem, Dist and
Rank, at word offsets 0 through 4 of its fixed frame.
Selection uses the fixed SELECTED slot for this cell.
A later cell's attempts cannot modify this frame.
The parent continuation of the first cell is
ANCHOR_EXTEND. All other parents are the CELL_NEXT
label of the immediately preceding cell. Targets here
are literal, so no dynamic stack pointer is required.

    CELL_INIT:
        u := extract cell i of Yj
        v := extract cell i of Pj
        f := extract cell i of C[j-1]
        key := add(PARITY_BASE[m],u)
        parity := load(key)
        s := xor(v,f)
        s := xor(s,parity)
        base := shl(s,m)
        base := add(LIST_BASE[m],base)
        store(Base,base)
        store(U,u)
        store(Rem,Remaining)
        store(Dist,Distance)
        store(Rank,0)
        fall through to CELL_NEXT

For the two systematic-only cells of field 1, replace CELL_INIT
above by extraction of the one systematic bit and the five stores,
using the literal SYS_ONLY_BASE as Base. The rest of this literal
block remains unchanged and its two entries correspond to e=0,1.

    CELL_NEXT:
        rank := load(Rank)
        stop := eq(rank,2^m)
        branch if stop to parent continuation
        base := load(Base)
        adr := add(base,rank)
        record := load(adr)
        cc := shr(record,17)
        needed := add(cc,62-DepthCell)
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
        d := shr(record,12)
        d := and(d,31)
        previous := load(Dist)
        Distance := add(previous,d)
        [at a field boundary, run FIELD or FINAL below]
        jump to next cell's CELL_INIT or,
            for the last cell, to its own CELL_NEXT

Here both possible exits make subsequent entries irrelevant:
there are no records left, or the sorted cost is already too
large. Loading a fresh rank and the parent path state at every
attempt restores that state after visiting any child subtree.
This is a depth-first traversal of exactly the admitted prefixes.
CELL_INIT costs at most 29 (18 to extract and 11 other
primitives). The single stopping attempt per initialized frame
costs at most 11. Each accepted attempt costs at most 28
including its final jump, exclusive of FIELD or FINAL.

At FIELD or FINAL, pack the 31 selected slices for this field
from their literal SELECTED addresses into Xj, using at most
31*(1+4)+6=161 operations. Compute Cj=Xj^3 using the cube routine.
At a nonfinal field, the next block can now use Cj. At FINAL
start TotalD from Distance and, for each of the six words of
AP and C2, load both, XOR, popcount and add. Compare TotalD
with BestD, then branch around replacement unless strictly
less. If replacing, load all 18 words of X and store them
in the corresponding B words, and move BestD to TotalD.
The anchor residual is P0+L(X0)+C2, so TotalD is exactly the
full distance: Distance already includes the anchor systematic
noise and both visited systematic/parity pairs. This
implements the specified distance ordering, with equality
retaining the first candidate. FIELD costs at most 180 plus
the cube and FINAL at most 250 plus the cube. All data live
in the indicated registers and arrays, including anchor Depth,
which is preserved through the child blocks. All shifts are
smaller than 64, all literal padding is zero, and all
operations on costs and distances fit in a word.

## Worst-case number of primitives

At cell index i, numbered from 1 through 62, a known parent remains
under budget r. This cell has a fixed width and an arbitrary syndrome.
For a systematic-only cell use the unique two-record row.
The only aspect of the row needed for the bound is the frequency
A_(s,c) of records of integer cost c, from the enumerated syndrome
weight histogram and the same cost function used by the fixed
decoder table. The checker considers all histograms that occur at
that width. At this parent the decoder accepts exactly the records
with c<=r-(62-i), by the stated monotone list traversal.

Define F_i(r) as a bound on the complete run of CELL_INIT and
all attempts at cell i, from the given r through return to its parent.
Set F_63(r)=0. Let T_i be zero for ordinary indices, CubeOps+180
for the first full-field boundary, and CubeOps+250 at i=62.
The exact recursive cost bound on a frame is

    F_i(r) = 40 + max_s sum_(c<=r-(62-i))
                         A_(s,c) * (28 + T_i + F_(i+1)(r-c)).

The 40 charge covers initialization and the unique stopping attempt.
Every accepted attempt charges 28, its field routine if any, and the
cost of the subsequent frame if any. Induct downward on i. At any
actual received word, field and choices so far, the frame's syndrome
is one of the enumerated rows. Downstream syndrome choices may be
fully adaptive, and each child budget is still exactly r-c, so
the independent maximum at its state remains valid for that child.
All recursive continuation jumps were included in the accepted attempt
bound or its stopping attempt. This induction proves the bound for
every input without requiring that a set of maximizing syndromes is
simultaneously realizable.

For any anchor depth d the frame starts at cap C(d).
Together with the increasing-subset traversal proof and the fixed
initial, parity setup, and field cube charges, an integer work bound is

    W = 1000 + sum_(d=0..R) binom(332,d) *
                 (44 + 410 + CubeOps + F_1(C(d))).

The verifier calculates all recurrence entries from 0 through
C(0) and checks the exact rational exponent inequality. The proof
uses the same caps for the word-RAM count and the error certificate.


The exact calculations validate the following rounded component
bounds and exponent with integer or rational operations.

    Unremoved full-code baseline error <= 55181057551535181/10^18
    Occupancy-weighted proper-run losses <= (373970, 395514)/10^18
    E_miss <= 1997643287089/10^18
    B = 27591527597795877/500000000000000000
    W = 1714363491736672601245918726983446372521051114591435656093119217826676119531910802664750335410883572974736894
    W^(500) <= 2^(179773).

Each length charge sums all starts. The error score adds the charges
and the miss bound to the unremoved full-code baseline.
