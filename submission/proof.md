# Five cyclic cube fields with Golay and Hamming parity

## Encoder and restricted list decoder

Use the ordered polynomial basis of
(K=\mathbb F_2[X]/(X^{200}+X^{30}+X^{19}+X^{11}+1)).
Write the 1000 message bits as five successive 200-bit field values
(b_0,\ldots,b_4), with the constant coefficient at the low end.
The codeword is the concatenation of these values followed by
five parity fields (p_i=L(b_i)+b_{i-1}^3) with cyclic indices.
Here (L) acts separately on the low-to-high sequence of sixteen
12-bit cells and two 4-bit cells. In a 12-bit cell form
(r=(uX^{11})\bmod(1+X^2+X^4+X^5+X^6+X^{10}+X^{11})).
Its low eleven parity bits are the coefficients of (r), and its high bit
is the parity of the sum of the Hamming weights of (u) and (r).
In a 4-bit cell use the linear map with columns
(14,13,11,7) as four-bit integers. These maps are linear. The
construction and all bit layouts are fixed in construction.py.
Irreducibility follows from the Rabin criterion for degree 200, whose
prime divisors are 2 and 5. The experiment applies repeated squaring
and checks the two gcds at indices 40 and 100 and the identity
(X^{2^{200}}=X\bmod P).

The decoder starts with the all-zero message as incumbent. It
enumerates all 200-bit configurations (x_0) at distance at most
(R) from the received systematic field (Y_0). For a given anchor,
process fields (j=1,2,3,4) in order. The candidate predecessor has
already been fixed. For each cell apply the following rule, with cell
width (m\in\{12,4\}), local parity map (L_m), systematic slice (u)
of (Y_j), and adjusted parity slice (v) of (P_j+x_{j-1}^3).
Put (s=v+L_m(u)). For every systematic error slice (e) put

    d_s(e) = wt(e) + wt(L_m(e)+s)
    A_s(e) = 9^(2m-d_s(e))
    S_s = sum_e A_s(e)
    c_s(e) = ceil(8 log2(S_s / A_s(e))).

The candidate cell is (u+e). The cost is evaluated in the fixed
table by bit_length of `(S_s^8-1)//A_s(e)^8`. Indeed for an integer
(n\ge0),
(S_s^8\le 2^n A_s(e)^8) exactly when
((S_s^8-1)\mathbin{//}A_s(e)^8<2^n).
All costs are positive because there is more than one positive term.
For each of the 72 cells the search traverses the entries in increasing
((d_s(e),e)), which gives nondecreasing cost. A new prefix of
length (\ell) is entered only if its total cost is at most
(C-(72-\ell)). A complete retained tuple has cost at most (C).
Conversely all prefixes of any tuple of cost at most (C) pass the
tests, since each future cost is at least 1. The final choice is the
first message of minimum total Hamming distance to the received word
among the initial incumbent and all complete tuples, using the traversal
order below. The anchor uses (R=33), precision is 8, and (C=1144).

## Missing the true word

For a true prefix all preceding candidate cubes equal the true cubes.
The true systematic and parity cell noises are therefore the original
independent channel errors. For a fixed syndrome (s) define
(h_s(d)=\#\{e:d_s(e)=d\}). For given (s), every (e) corresponds to
exactly one pair of noise slices, with weight (d_s(e)).
The probability-generating polynomial of the true cost in a cell is

    (sum over s,d of h_s(d)*9^(2m-d)*Z^c_s(d)) / 10^(2m).

The true anchor systematic noise depends on disjoint channel bits
from all 72 child cells. The product of 64 factors of width 12 and
8 factors of width 4, kept through degree (C), gives the chance
that the true child costs sum to at most (C). Multiply this by the
exact 200-bit binomial probability through (R). Its complement
bounds the probability (E_{\rm miss}) that the true word is absent.
A true codeword that equals the fallback only improves this bound.
The experiment calculates the cell histograms by complete enumeration
and groups equal histograms with their integer multiplicities.
It multiplies truncated polynomials by packing each coefficient
in a 5400-bit digit. Every retained coefficient at every step is
bounded by (10^{1600}<2^{5400}), the total untruncated channel
numerator on the 1600 child bits, so no carry affects an adjacent
coefficient.

For a fixed received word, take full-code maximum likelihood (ML)
with a tie preference for the restricted decoder's output when
that output is among the full-code minimizers. ML success probability
under the uniform message is the same for any tie preference because
all minimizers have equal likelihood at that received word.
Whenever this comparator is correct and the true word is in the
restricted set, the restricted decoder sees a full-code minimizer
and must return one. The comparator then uses that restricted
output and is correct only if the restricted output is true.
Consequently restricted block error is at most the full-code ML
block error plus (E_{\rm miss}).

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
difference in field i is L(e_i) plus the cube
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
The checker verifies D_t is nondecreasing through t=230.

## Competitors changing a proper subset

Let (q=2^{200}) and (Q=q^5). For a proper nonempty set of changed
systematic fields, its cyclic runs of consecutive changes have disjoint
affected coordinates. A run (i,\ldots,i+a-1) affects its (a)
systematic fields and the (a+1) parity fields (i,\ldots,i+a).
The distance change of the whole competitor is the sum of the changes
on these disjoint sets. Thus any proper competitor of distance at
most that of a given true word implies a single-run competitor of
distance at most that of that word. We bound removal for inclusive
comparisons, including all ties.

Fix the true noise and a run of length (a\in\{1,2,3,4\}).
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
jointly have (u=200(2a-1)) bits. Each potential systematic error
pattern fixes the direction on the run. We can upper bound the
expected number of competitors of local distance at most (w) by
((2/q)^a) times the coefficient prefix through (w) of

    (1+Z)^u * product_(18 pair cells) H_s(Z).

Here (w) is the local true noise weight in all (400+u) affected
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

    2^(-floor(199*a/2)) * r^(-w) * (1+r^2)^(u/2)
        * product_(18 pair cells) sqrt(H_s(r^2)).

This rounds the power of two upward. Use rational
(r=46/100,44/100,42/100,40/100) for lengths 1, 2, 3, 4 respectively.
The verifier bounds each square root from above at resolution
(2^{-64}). Write an evaluated (H_s(r^2)=n/d) and take (j)
to be the integer square root of the floor of (2^{128}n/d).
Then the exact numerator of the rounded bound is (j) if
(j^2d=2^{128}n), or (j+1) otherwise.

For each length, to average over the channel noises, form the
polynomial of each paired cell with coefficient at degree (d)
equal to the sum of (h_s(d)) times the rounded integer root
numerator over all syndromes (s). Multiply the sixteen width-12
and two width-4 polynomials and ((1+Z)^u). At local degree (w)
apply (r^{-w}), via a common denominator using the largest local
degree that is needed, 230. Multiply by ((1+Z)^{2000-400-u})
to include outside weights. A final coefficient at total degree
(v) is multiplied by the common run factors and the channel
factor (9^{2000-v}/10^{2000}). Denote this bound on the joint
run-removal probability at weight (v) by (J_{a,v}). The
translation-invariant channel makes the bound identical at each of
the five starts. Terms above total degree 230 are never needed.

In this computation the rounded root numerator is at most (2^{71}):
(H_s(r^2)\le H_s(1)=2^m) and (m\le12).
The sum of the integer coefficients of the paired factors
is at most (2^{400+18\cdot71}), and the binomial factor contributes
(2^u\le2^{1400}). Thus the first product is safe in 4000-bit digits.
Writing (r=a_0/b_0), the local weight multiplier on the integer
coefficients is (b_0^w a_0^{230-w}\le100^{230}<2^{1529}).
After the remaining binomial factor the resulting sum is bounded by
(2^{2000+18\cdot71+1529}<2^{6000}). Truncations and binary powers
therefore do not create carries in the exact product.

## Full-code occupancy

At distance cutoff (t\le230), remove the pairs of true word and
noise whose noise weight is at most (t) and for which a single-run
competitor has distance no greater. Write (L_t) for the number of
these pairs over all (Q) messages. Two retained candidates at a
common received word must then differ in all five fields: otherwise
their proper-set comparison would make the farther one removable
(or both if equally far). There are (Q V_t-L_t) remaining
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
            + 5*sum_(a=1..4) sum_(v<=230) alpha_v*J_(a,v),

with (U_{-1}=0). The omitted contribution to success beyond
the cutoff is nonnegative. Both bounds used in the final error
score are checked with integer and rational arithmetic.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The ten received fields are copied to padded little-endian
four-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, five frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed tables then contain reduction, square, local parity,
and sorted candidate rows. The Golay table has 2^24 words, and
the Hamming table 2^8 words. Each word is an error in bits 0..11,
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
can be discarded. All 72 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use four accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

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
Use this for the 40 padded words of Y and P, then store each.
Load, popcount and add each padded word to form the distance
of the all-zero incumbent, at cost 120. Store zero in the 20
padded incumbent words. Copy Y0 to X0 in eight primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear four
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the four
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=200, and Depth is a persistent register.

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
N=sum_(i<=R) binom(200,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base, load the parity image, and pack
these images into a four-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Move Remaining=C and
Distance=Depth, and jump to the first CELL_INIT below.
This uses one cube and at most 250 further primitives:
the 18 cells use at most 18*(6+2+4), and finishing AP
uses 12 plus the single jump.

### Cell traversal

Emit one literal block for each of the 72 cells, in
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

    CELL_NEXT:
        rank := load(Rank)
        stop := eq(rank,2^m)
        branch if stop to parent continuation
        base := load(Base)
        adr := add(base,rank)
        record := load(adr)
        cc := shr(record,17)
        needed := add(cc,72-DepthCell)
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

At FIELD or FINAL, pack the 18 selected slices for this field
from their literal SELECTED addresses into Xj, using at most
18*(1+4)+4=94 operations. Compute Cj=Xj^3 using the cube routine.
At a nonfinal field, the next block can now use Cj. At FINAL
start TotalD from Distance and, for each of the four words of
AP and C4, load both, XOR, popcount and add. Compare TotalD
with BestD, then branch around replacement unless strictly
less. If replacing, load all 20 words of X and store them
in the corresponding B words, and move BestD to TotalD.
The anchor residual is P0+L(X0)+C4, so TotalD is exactly the
full distance: Distance already includes the anchor systematic
noise and all four visited systematic/parity pairs. This
implements the specified distance ordering, with equality
retaining the first candidate. FIELD costs at most 180 plus
the cube and FINAL at most 200 plus the cube. All data live
in the indicated registers and arrays, including anchor Depth,
which is preserved through the child blocks. All shifts are
smaller than 64, all literal padding is zero, and all
operations on costs and distances fit in a word.

## Worst-case number of primitives

Fix a received word and a visited anchor. In the unrestricted
tree of child choices, assign each edge at its prefix the
weight A_s(e)/S_s. These weights sum to 1 for every parent,
including when its syndrome depends on previously selected
fields. The product weights over all prefixes of a fixed
length sum to 1 by induction. A visited prefix of length
ell has product weight at least
2^(-(C-(72-ell))/8) by the definition of its cost.
Thus the number of visited prefixes at that level is at most
n_ell=2^ceil((C-(72-ell))/8), and put n_0=1.
For level ell there are at most n_(ell-1) initializations
and failed attempts, and n_ell successful attempts. When ell
is a positive multiple of 18, there are at most n_ell field
boundary routines, including scoring at ell=72.

Let CubeOps=S+Mmul from the arithmetic macros, let F_ell
be zero at nonmultiples of 18, CubeOps+180 at the first three
field boundaries, and CubeOps+200 at the last. A bound
valid on every received word is

    W = 1000 + N * (44+250+CubeOps
          + sum_(ell=1..72) (40*n_(ell-1) + (28+F_ell)*n_ell)).

The computation verifies the declared exact exponent against
this integer. The proof uses the same capped search in the
admission argument and in the word-RAM count.

For this profile the exact computed bounds are

    E_ML <= 56052834598407800/10^18
    E_miss <= 4833234674510508/10^18
    B = 15221517318229577/250000000000000000
    W = 3621989781858856031893530258219563868744481202792556581271415974234845706154037234896
    W^(1000) <= 2^(280899).
