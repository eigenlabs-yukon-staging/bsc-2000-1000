# Ten cyclic cube fields with Golay and Hamming parity

## Encoder and restricted list decoder

Use the ordered polynomial basis of
(K=\mathbb F_2[X]/(X^{100}+X^{29}+X^{23}+X^{16}+1)).
Write the 1000 message bits as ten successive 100-bit field values
(b_0,\ldots,b_9), with the constant coefficient at the low end.
The codeword is the concatenation of these values followed by
ten parity fields (p_i=L(b_i)+b_{i-1}^3) with cyclic indices.
Here (L) acts separately on the low-to-high sequence of eight 12-bit cells and
one 4-bit cell. In a 12-bit cell form
(r=(uX^{11})\bmod(1+X^2+X^4+X^5+X^6+X^{10}+X^{11})).
Its low eleven parity bits are the coefficients of (r), and its high bit
is the parity of the sum of the Hamming weights of (u) and (r).
In a 4-bit cell use the linear map with columns
(14,13,11,7) as four-bit integers. Both maps are linear. The
construction and all bit layouts are fixed in construction.py.
Irreducibility follows from the Rabin criterion for degree 100, whose prime
divisors are 2 and 5. The experiment applies repeated squaring
and checks the gcds at indices 20 and 50 and the identity
(X^{2^{100}}=X\bmod P).

The decoder starts with the all-zero message as incumbent. It
enumerates all 100-bit configurations (x_0) at distance at most
(R) from the received systematic field (Y_0). For a given anchor,
process fields j=1 through 9 in order. The candidate predecessor has
already been fixed. For each cell apply the following rule, with cell
width (m\in\{12,4\}), local parity map (L_m), systematic slice (u)
of (Y_j), and adjusted parity slice (v) of (P_j+x_{j-1}^3).
Put (s=v+L_m(u)). For every systematic error slice (e) put

    d_s(e) = wt(e) + wt(L_m(e)+s)
    A_s(e) = 9^(2m-d_s(e))
    S_s = sum_e A_s(e)
    c_s(e) = ceil(16 log2(S_s / A_s(e))).

The candidate cell is (u+e). The cost is evaluated in the fixed
table by bit_length of `(S_s^16-1)//A_s(e)^16`. Indeed for an integer
(n\ge0),
(S_s^16\le 2^n A_s(e)^16) exactly when
((S_s^16-1)\mathbin{//}A_s(e)^16<2^n).
All costs are positive because there is more than one positive term.
For each of the 81 cells the search traverses the entries in increasing
((d_s(e),e)), which gives nondecreasing cost. A new prefix of
length (\ell) is entered only if its total cost is at most
C(d)-(81-\ell). A complete tuple has cost at most C(d).
Conversely all prefixes of any tuple of cost at most C(d) pass the
tests, since each future cost is at least 1. The final choice is the
first message of minimum total Hamming distance to the received word
among the initial incumbent and all complete tuples, using the traversal
order below. The anchor radius is R=23 and cost precision is 16. For each
anchor distance d from 0 through R use C(d)=2715-48d.

## Missing the true word

For a true prefix all preceding candidate cubes equal the true cubes.
The true systematic and parity cell noises are therefore the original
independent channel errors. For a fixed syndrome (s) define
(h_s(d)=\#\{e:d_s(e)=d\}). For given (s), every (e) corresponds to
exactly one pair of noise slices, with weight (d_s(e)).
The probability-generating polynomial of the true cost in a cell is

    (sum over s,d of h_s(d)*9^(2m-d)*Z^c_s(d)) / 10^(2m).

The true anchor systematic noise depends on disjoint channel bits
from all 81 child cells. The product of 72 width-12 and 9 width-4
factors has the distribution of their total true cost.
For each anchor distance d from 0 through R, multiply its exact
binomial mass binom(100,d)*9^(100-d)/10^100 by the product's
coefficient prefix through C(d). The sum is the actual true-path
admission probability before including the fallback. Its complement
bounds the error charge for missing the true path.
A true codeword that equals the fallback only improves this bound.
The experiment calculates the cell histograms by complete enumeration
and groups equal histograms with their integer multiplicities.
It multiplies truncated polynomials by packing each coefficient
in a 6000-bit digit. Every retained coefficient at every step is
bounded by (10^{1800}<2^{6000}), the total untruncated channel
numerator on the 1800 child bits, so no carry affects an adjacent
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

## Competitors changing all ten fields

Write q=2^100 and Q=q^10. For nonzero e, the cube
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
annihilate trace on a spanning set. The cubic root bound gives at most three nonzero directions e
for each nonzero normal mask.

A full-change competitor is specified by nonzero
directions e_i in all ten fields. Its parity
difference in field i is L(e_i) plus the cube
difference in field i-1. Thus, for fixed directions,
the ten translated hyperplane distributions are
independent over the uniform true information fields.
The probability of a radius-l ball in the 1000
parity bits can be evaluated by expanding the
ten affine indicators. Multiplication by their
common point mass 1024/Q leaves 1/Q times a sum
of 1024 terms. The empty mask term is the
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

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(100*h,w).

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
most binom(100*k,a) choices of the relevant
directions. With total systematic error weight b
on the other fields there are at most
binom(100*(10-k),b) outside choices. Thus for a
candidate total distance threshold t, the sum
of this term over directions with weights a,b,
where a+b<=t, is bounded by

    binom(100*(10-k),b) *
    H_k(t-a-b, binom(100*k,a)).

This accounts for any fixed actual noise because
a candidate error on a systematic field determines
its direction uniquely. Sum these bounds over a,b
and the binom(10,k) selections of masks, and add
the 2000-bit ball volume for the empty term. The
empty term is the sum of the parity ball volume
weighted by all systematic binomial counts.
Call the resulting integer D_t. The expected
number of full-change competitors within total
t for each fixed true noise is at most D_t/Q.
The checker verifies D_t is nondecreasing through t=230.

## Competitors changing a proper subset

Let q=2^100 and Q=q^10. For removal, order candidates at each
received word first by distance and then by their numerical systematic
message. A proper nonempty set of changed systematic fields breaks
into cyclic runs. One run of length a from i through i+a-1 affects its
a systematic fields and the a+1 parity fields from i through i+a.
Distinct runs affect disjoint coordinates. If the whole competitor
is closer, some run decreases distance. If the whole competitor
ties and no individual run decreases distance, all run changes are
zero. The run containing the highest changed message bit is ordered
first by the numerical tie rule exactly when the whole competitor
is ordered first. Hence any preferred proper competitor implies a
preferred single-run competitor.

Fix a run of length a and fix the true noise. Its first systematic
and parity pair has changes e and L(e), since the predecessor is
unchanged. For each cell of that pair the exact syndrome is invariant.
Let H(Z) be the product of these nine fixed syndrome polynomials
H_s(Z)=sum_d h_s(d) Z^d. The possible errors of the first pair with
given syndromes are counted by H, and their systematic values identify
the first direction uniquely. Outside this pair, the affected
coordinates contain n=(a-1)*100 systematic bits and p=a*100 parity
bits. Put u=n+p. For an actual fixed run direction the a parity
differences outside the pair are independent translated affine
hyperplanes with mass 2/q per available word, using the distinct
true predecessor fields. In particular the expected number of strict
and inclusive competitors within any local total threshold is at
most (2/q)^a times the corresponding strict and inclusive prefixes
of H(Z)*(1+Z)^u. Enlarging the count by allowing extra systematic
patterns does not assign hyperplane mass to zero directions, since
each counted fixed actual direction is nonzero before enlargement.

### Ranked bound outside the first pair

For any fixed nonzero first direction and fixed true noise, we also
bound the expected number of settings of the other directions whose
total candidate error outside the first pair is at most z. Write
D_(a,z)/q^a for this bound. Expand the parity-ball membership of
their a independent translated hyperplanes, as in the full-change
analysis but with p parity bits and n systematic bits. The empty
mask contribution to D_(a,z) is the ordinary Hamming ball volume
on p+n bits. For a nonempty selection of masked predecessor fields
let k be the number among the other a-1 predecessors. There are
binom(a-1,k) such selections for each choice of whether the fixed
first mask is also included. Write K_l(w) for the signed p-bit
parity ball sum at radius l and combined mask weight w. The earlier
Krawtchouk recurrence holds with p in place of 1000.

For each k>0 the other selected masks have weight-w capacities

    3^k * sum_(h=1..k) (-1)^(k-h)*binom(k,h)*binom(100*h,w).

If the first mask is excluded, assign these entries magnitude
abs(K_l(w)). If it is included, assign them magnitude
max_(1<=j<=100) abs(K_l(w+j)).
The latter is uniform over the first direction. At a candidate
systematic weight h on the k selected predecessor fields, at
most binom(100*k,h) choices of their directions are possible.
Use the sum of that many largest magnitudes at the listed
capacities, just as for the full-change envelope. The unselected
other systematic fields contribute at most binom(n-100*k,b)
choices at their candidate weight b. Sum over h,b with h+b<=z at
l=z-h-b and multiply by binom(a-1,k) for each k and mask choice.
For k=0 the first-mask-only selection uses h=0 and
max_(1<=j<=100) abs(K_l(j)); the empty selection is already covered.
This defines D_(a,z) independently of the first direction and all
true noise. The verifier computes the sliding maxima of the signed
rows for the additional first mask and the grouped top sums for
the other masks. There is no dependence on the mask's translation
sign, because every nonempty term uses an absolute upper bound.

For a fixed first pair error pattern of weight d its outside
threshold is z=w-d, where w is the whole local true noise weight.
Sum D_(a,w-d)/q^a over the possible first pair errors. If the first
direction for a pattern is zero, omit it in the true count and
still allow that nonnegative expression in the upper bound. There
is consequently a strict competitor expectation bound using
[Z^(w-1)] H(Z)*D_a(Z)/q^a and an inclusive bound using
[Z^w] H(Z)*D_a(Z)/q^a, where D_a(Z) has coefficient D_(a,z)
at z from 0 to 230, and coefficients for negative z vanish.

### Tie averaging, clipping and averaging true weights

Fix exact first-pair syndromes, the local true noise weight w
and the complete outside noise vector. All true words and noise
vectors of that class have equal weight under the uniform message
and BSC. Every equal-distance directed pair has a reverse directed
pair at the same received word. This involution preserves the exact
syndromes, the local weight and every outside bit. The two distinct
messages have exactly one preferred numerical ordering, so the
class-average count of preferred comparisons is half the sum of the
strict and inclusive counts. Clipping this average count above by
one bounds the average indicator that any preferred competitor exists.

Let A_w=[Z^w]H(Z)*(1+Z)^u, V_w=sum_(j<=w) A_j,
and let T_w=[Z^w]H(Z)*(1+Z)*D_a(Z). The two methods above
therefore give the class-average probability upper bound

    min(2*q^a, 2^a*(2*V_w - A_w), T_w) / (2*q^a).

There are A_w local true noise patterns in that class, so the
term is multiplied by A_w and summed over the exact syndromes.
Because the expression depends only on their profiles, group equal
Golay and Hamming cell polynomials by multiplicities. For the eight
Golay positions use the multinomial coefficient for a collection
of profile types and the product of their individual syndrome
multiplicities. Then include the final Hamming cell with its
profile multiplicity. Finally multiply by the outside binomial
coefficient at weight v-w and the channel factor
9^(2000-v)/10^2000, and sum over w to bound J_(a,v), the
run-removal probability at total weight v. The result applies to
each of the ten starts for each length a=1 through 9.

Each H_s has only one degree parity, as the checker verifies from
its local enumerator. For packed products, shift each H_s to its
first nonzero degree and use its half exponent variable, then
multiply by the even and odd parts of the unpaired polynomial,
and independently by those of (1+Z)*D_a(Z). Each coefficient
used through degree 230 is bounded by 2^100 times the larger
of the sums of the two full input coefficient vectors through
degree 230. This follows because H(1)=2^100 for each choice of
syndromes. Therefore the bit length of that product, rounded up
to full bytes, suffices for every packed digit. All products are
nonnegative, so terms in higher digits cannot affect the extracted
lower ones. The declared experiments check these exact integers
for every run length.

## Full-code occupancy

At distance cutoff (t\le230), remove the pairs of true word and
noise whose noise weight is at most (t) and for which a preferred
single-run competitor exists. Write (L_t) for the number of
these pairs over all (Q) messages. Two retained candidates at a
common received word must then differ in all ten fields: otherwise
their proper-set comparison would make the one later in the removal order removable. There are (Q V_t-L_t) remaining
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
            + 10*sum_(a=1..9) sum_(v<=230) alpha_v*J_(a,v),

with (U_{-1}=0). The omitted contribution to success beyond
the cutoff is nonnegative. Both bounds used in the final error
score are checked with integer and rational arithmetic.

## Word-RAM program

The counted program uses the exact addresses in construction.py.
The twenty received fields are copied to padded little-endian
two-word arrays Y and P. Candidate fields X, their cubes C, and
the incumbent B use the same format. AP holds the anchor parity
residual without the last field cube. Other working regions store
the anchor traversal rows NEXT and LAST, one SELECTED word per
child cell, five frame words per child cell, and cube scratch.
The arrays are consecutive after the prescribed input/output.
The fixed tables then contain reduction, square, local parity,
and sorted candidate rows. The Golay table has 2^24 words
and the Hamming table has 2^8 words. Each word is an error in bits 0..11,
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
can be discarded. All 81 DFS frames use RAM, and the cell and
anchor routines each need fewer than 20 scratch registers.
The cube routines use two accumulators and a fixed number of
temporary registers below 20. These scratch slots are reused,
so no path needs more than 64 registers.

Use k=2 words and h=36 high-word bits for a field.
For literal s=1 and 4, the reduction table REDUCE_s
has 2^s entries, with entry c the carryless product
of c and P+X^100. This fits in one word since the
largest possible degree is 32. To multiply a field
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
c*X^100 by its residue gives the last XOR. This
SHIFT_s macro costs 3k+3 primitives.

To square a literal source array, start k zero
accumulators. For each of b0=13 byte slots, emit a
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
Move R registers to zero. For each of n0=25 nibble
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
padded incumbent words. Copy Y0 to X0 in four primitives,
set Depth to zero and store NEXT[0]=0.

For packing a field from known short pieces, clear two
registers by free moves and insert a piece with a left shift
and OR. On a crossing also right shift the original piece
and OR it into the next word. No dynamic tests occur because
every index and width is literal. At the end store the two
accumulators. For output, each literal output word draws
pieces from at most three successive padded incumbent words.
Even extracting these one by one and inserting them uses
at most 3*(6+2)+1 primitives per word. The final partial
output word is padded by the initialized accumulator. The
sum of startup and output uses fewer than 1000 operations.

### Anchor tree

The anchor starts at X0=Y0. Its increasing subsets of
flipped positions use the following template. R is the
literal radius, M=100, and Depth is a persistent register.

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
N=sum_(i<=R) binom(100,i) visited anchors. These events
cost at most 23, 6, 15 and 2 respectively, so their total
cost is at most 44N.

At ANCHOR_VISIT compute C0=X0^3 using the cube routine.
To prepare AP, extract each literal cell of X0, add its
local parity table base, load the parity image, and pack
these images into a two-register field accumulator as
above. For each accumulator load the corresponding P0
word, XOR and store to AP. Set Remaining=C(Depth) with one multiplication by 48 and one
subtraction from 2715. Move Distance=Depth and jump to CELL_INIT.
This uses one cube and at most 250 further primitives:
the 9 cells use at most 9*(6+2+4), and finishing AP
uses 6, two budget operations, and one jump.

### Cell traversal

Emit one literal block for each of the 81 cells, in
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
        needed := add(cc,81-DepthCell)
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

At FIELD or FINAL, pack the 9 selected slices for this field
from their literal SELECTED addresses into Xj, using at most
9*(1+4)+2=47 operations. Compute Cj=Xj^3 using the cube routine.
At a nonfinal field, the next block can now use Cj. At FINAL
start TotalD from Distance and, for each of the two words of
AP and C9, load both, XOR, popcount and add. Compare TotalD
with BestD, then branch around replacement unless strictly
less. If replacing, load all 20 words of X and store them
in the corresponding B words, and move BestD to TotalD.
The anchor residual is P0+L(X0)+C9, so TotalD is exactly the
full distance: Distance already includes the anchor systematic
noise and all nine visited systematic/parity pairs. This
implements the specified distance ordering, with equality
retaining the first candidate. FIELD costs at most 180 plus
the cube and FINAL at most 200 plus the cube. All data live
in the indicated registers and arrays, including anchor Depth,
which is preserved through the child blocks. All shifts are
smaller than 64, all literal padding is zero, and all
operations on costs and distances fit in a word.

## Worst-case number of primitives

At cell index i, numbered from 1 through 81, a known parent remains
under budget r. This cell has a fixed width and an arbitrary syndrome.
The only aspect of its row needed for the bound is the frequency
A_(s,c) of records of integer cost c, from the enumerated syndrome
weight histogram and the same cost function used by the fixed
decoder table. The checker considers all histograms that occur at
that width. At this parent the decoder accepts exactly the records
with c<=r-(81-i), by the stated monotone list traversal.

Define F_i(r) as a bound on the complete run of CELL_INIT and
all attempts at cell i, from the given r through return to its parent.
Set F_82(r)=0. Let T_i be zero for ordinary indices, CubeOps+180
for the first eight full-field boundaries, and CubeOps+200 at i=81.
The exact recursive cost bound on a frame is

    F_i(r) = 40 + max_s sum_(c<=r-(81-i))
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

    W = 1000 + sum_(d=0..R) binom(100,d) *
                 (44 + 250 + CubeOps + F_1(C(d))).

The verifier calculates all recurrence entries from 0 through
C(0) and checks the exact rational exponent inequality. The proof
uses the same caps for the word-RAM count and the error certificate.

The four declared runs check the following integer data for the
common denominator 10^18 and the integer work inequality.

    Unremoved full-code baseline error <= 58309282338755975/10^18
    Occupancy-weighted proper-run losses <= (32130851529865, 654771251756452, 1530340318922916, 1915670610504305, 1775789974200312, 1300780205703891, 705126443021130, 210429773308846, 8842704647383)/10^18
    E_miss <= 19972840061346590/10^18
    B = 31924188747210713/200000000000000000
    W = 70484566301335116987477129842234353111923230085460388110258
    W^(100) <= 2^(19549).

The claimed error adds ten times the sum of run charges and
the miss charge to the unremoved full-code baseline.
