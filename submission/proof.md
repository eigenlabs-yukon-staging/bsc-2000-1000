# Six search charts for a two-coordinate field code

## Channel code and actual decoder

Let m=500, K=F2[X]/(P) with P=X^500+X^21+X^13+X^8+1,
and q=2^m. Encode field elements in the polynomial basis, with
coefficient of X^i in bit i. Write Q=q^2. Split the 1000-bit
message into elements a,b in that order, and encode it in four
successive 500-bit groups by

    f(a,b) = (a,b,c,d),   c=a*b,   d=a^3+lambda*b^3,

where lambda is the field element represented by integer 3.
This encoder is defined for all messages and is injective
because it has both message groups as systematic coordinates.

The declared computation checks X^(2^500)=X modulo P and
gcd one with both X^(2^250)+X and X^(2^100)+X by repeated
modular squaring and exact polynomial division. Indeed an
irreducible factor of degree j divides X^(2^l)+X exactly
when j divides l. In its quotient field the j-th Frobenius
power fixes all elements by the order of the multiplicative
group. No smaller positive power fixes all elements, by
the root bound for a polynomial of degree less than 2^j.
Fixing the residue of X fixes the whole quotient since it
generates that field. Thus the first check forces all
factor degrees to divide 500. Every proper divisor of 500
divides 250 or 100, so the gcd checks exclude all proper
degrees. This proves irreducibility. Set h=(q-1)/3 and
zeta=lambda^h. The experiment checks that zeta has order
three. In particular lambda is not a cube, because a cube
x^3 has (x^3)^h=1 for x!=0.

The decoder considers only messages with a,b both nonzero.
An eligible message must also have some pair of its four
encoded groups at total distance at most R=108 from the
corresponding received groups. Among all eligible messages
it chooses one minimizing total distance across the four
groups, using the smallest integer a+2^500 b on ties.
The set is always nonempty. On the systematic pair, any
zero received group can be replaced by field element 1,
using at most two flips in total. Every other received
systematic group can be copied. The decoder's exclusion
of zero field components is part of its definition and
will be reflected in the success analysis.

## Parity differentials and multiplicities

For a fixed nonzero systematic difference (s,t), the first
parity difference is a*t+b*s+s*t. The second is

    s*a^2+s^2*a+s^3+lambda*(t*b^2+t^2*b+t^3).

These are affine F2-linear functions of a,b. Suppose s!=0.
Fixing the first parity difference leaves b affine in a with
its linear coefficient t/s. Substituting into the second
gives an affine function of a whose linear part is
A*(a^2+s*a), where A=(s^3+lambda*t^3)/s^2.
This scalar is nonzero, since if t!=0 an equality to zero
would make lambda a cube, and if t=0 it equals s.
The kernel is {0,s}. If s=0 then t!=0. Fixing the first
parity difference determines a, while b maps into the
second via lambda*(t*b^2+t^2*b) with kernel {0,t}.
So in every direction each occupied parity difference
has exactly two preimages. The parity difference for
a uniform message is uniform on an affine hyperplane
in the 1000-dimensional binary parity space.

We bound how often different directions can have the same
binary mask for this hyperplane. The F2 trace on K is
Tr(x)=sum_(i=0..499) x^(2^i). It is F2-linear with values
in F2, and nonzero by the root bound. Its coordinate pairing
is nonsingular. If Tr(z*x) were zero for all x for z!=0,
the trace would vanish everywhere since multiplication by z
is surjective. Thus every binary parity mask has a unique
trace pair (alpha,beta) such that its form on (c,d) is
Tr(alpha*c+beta*d).

The mask annihilates the linear part of the derivative in
direction (s,t) precisely when the two coefficients of
a and b, after absorbing squares under the trace, vanish:
    
    alpha*t + beta*s^2 + sqrt(beta*s) = 0,
    alpha*s + beta*lambda*t^2 + sqrt(beta*lambda*t) = 0.

Here the field square root is unique, and the identity
Tr(z*a^2)=Tr(sqrt(z)*a) follows by squaring under the trace.
If beta=0 and the mask is nonzero, these equations force the
zero direction. If alpha=0 and beta!=0, squaring the
equations gives one nonzero quartic each in s and t,
so there are at most sixteen pairs, including zero.
When alpha and beta are nonzero, square the equations.
The first expresses t^2 as A0*s^4+B0*s with A0!=0.
Substituting its square into the second expresses t as
A1*s^8+B1*s^2 with A1!=0. Squaring and comparing with
the first expression for t^2 gives a nonzero degree
sixteen polynomial in s, and t is uniquely fixed by s.
Thus the total solution count is at most sixteen including
zero. A given nonzero mask can therefore occur for at
most fifteen nonzero directions. The hyperplane for a
given direction has a unique nonzero mask because of its
proved rank.

## Exact coverage bound

For each total error weight t, let E_t consist of the errors
of weight at most t with some pair of groups of combined
weight at most R. Let V_t=|E_t|. For a received vector y let
m_y count good messages, meaning both field components
nonzero, with y+f(a,b) in E_t. Write G=(q-1)^2 for the count
of such messages. Then sum_y m_y=G V_t.

For a fixed true error z in E_t consider the number of
competing messages, dropping the good-message restriction
when bounding pairs. For a fixed systematic difference
(s,t') its two systematic error groups are the two fields
of that difference XOR the corresponding groups of z.
This is independent of the starting message. Denote the
two systematic error weights by u,v and put ell=t-u-v.
If u+v<=R, the permitted parity error weights p_0,p_1
satisfy simply p_0+p_1<=ell. Otherwise the condition is

    p_0+p_1<=ell and
    (p_0+p_1<=R or min(p_0,p_1)<=R-min(u,v)).

For the latter, a mixed pair is possible exactly if the
smaller parity and systematic weights sum to at most R.
The parity sets in both cases are symmetric under swapping
the two parity groups and invariant under permutations
within each. Their volume J is computed by summing
binom(m,p_0)binom(m,p_1) over the specified weight pairs.

Let D=42. A half-mask of weight between D and m-D inclusive
is balanced. For that half, its signed Hamming-shell
coefficient at weight w has absolute value at most
C_w=[x^w](1+x^2)^D(1+x)^(m-2D). In fact, its exact
signed sum is the coefficient of (1-x)^s(1+x)^(m-s) for
mask weight s. Pair at least D opposite factors and take
absolute values coefficientwise. Surplus pairs are bounded
by (1+x)^2 in place of 1+x^2. For the other half the
trivial bound is its shell cardinality. If either half
is balanced, the absolute signed sum over the parity
set is at most

    C = sum over its permitted p_0,p_1 of C_(p_0)*binom(m,p_1).

If the second half is the balanced one, swap the parity
weights in this sum using the set's symmetry. Translating
the set by the fixed actual parity noise multiplies the
signed sum by a sign. An affine hyperplane has either
sign in the intersection formula, so its intersection
with the translated parity set is at most (J+C)/2.
If both halves are unbalanced the trivial bound is J.
There are at most
Kbad=15*(2*sum_(w=0..D-1) binom(m,w))^2
such directions. This counts both extreme ranges by
complementation and uses the proved multiplicity limit.
The coefficient bound also ensures 0<=C<=J.

Group the systematic weight pairs by total u+v and their
allowed parity shape. The shape is either the unrestricted
ball up to ell or the one determined by R-min(u,v)
clamped below by -1. For a class let N be the sum of
binom(m,u)binom(m,v) over its ordered pairs. This upper
bounds its number of nonzero systematic differences,
since translating by the actual systematic error is a
bijection. The expected competitor count in the class
over all uniform starting messages is at most

    [N*(J+C)+min(N,Kbad)*(J-C)] / Q.

This uses the mass 2/Q per value on the affine hyperplane.
Any good-message restrictions can only decrease the
ordered distinct-pair count. Sum the numerators over
these classes to get A_t. The exact volume is V_t=sum N*J
over the same classes, since in counting all admissible
errors the two systematic weights are unrestricted except
for the displayed shape. The declared computation implements
these sums. It precomputes prefix sums on parity shells.
At parity total w<=R all the shell is included. When
w>R it includes the union p_0<=R-min(u,v) or
p_0>=w-(R-min(u,v)); its prefix/suffix sum counts an
overlap only once by taking the whole shell when the two
ranges meet. This works both for J and for the coefficient
bound C.

Sum over true errors z in E_t and all uniform starting
messages to bound the total ordered pair count by V_t A_t.
Indeed, each ordered pair of different candidates at a
received y corresponds to one true message, that error z,
and its competitor at the same y. Counting over all maps
with no restriction on good messages in this bound can only
increase the count for m_y. Thus
sum_y m_y^2<=V_t*(G+A_t). Set F_t equal to the number of y
with positive m_y divided by Q. Cauchy-Schwarz gives
F_t>=V_t G^2/(Q(G+A_t)).
For a positive integer occupancy x and integer k>=1,
(x-k)(x-k-1)>=0 gives k(k+1)>=(2k+1)x-x^2.
This extends to the zero occupancy with zero left side.
Summing yields
F_t>=V_t*(2kG-A_t)/(Q*k*(k+1)).
The computation uses k=1+floor(A_t/G), floors the two bounds
and keeps their maximum with the previous lower bound.
Write its nondecreasing sequence as D_t.

For a received word with winning distance w, its contribution
to joint average success, summed over all Q uniform messages,
is Q^-1 * 9^(2000-w)/10^2000. There is exactly one selected
message; its noise vector is determined and ties have the
same likelihood. The count at winning distance w is exactly
Q*(F_w-F_(w-1)). Consequently summation by parts gives the
lower success numerator
sum_(t=0..240) (D_t-D_(t-1))*9^(2000-t).
The coefficient of each cumulative D_t is nonnegative,
so replacing those numbers by the actual F_t can only
increase the sum. This counts exact success through t=240,
with further weights nonnegative. The run checks that one
minus this numerator over 10^2000 is at most
27654689647274881/500000000000000000. Normalization is
by all Q messages, including the ones the decoder excludes.

## The six field charts

The chart for a specified pair of group values produces
exactly its compatible (a,b) with a,b nonzero. Chart pairs
are visited in lexicographic pair order. On each chart
enumerate all two-group received perturbations of weight at
most 108. Thus the union of the output messages is exactly
the eligible set in the decoder definition. Duplicate messages
do not alter the minimum or numerical tie rule.

Pairs (a,b), (a,c) and (b,c) determine the message directly:
skip keys with a required nonzero value zero, and for
(a,c) or (b,c) recover the other component by division.
For (a,d) find all roots of b^3=lambda^-1*(d+a^3).
For (b,d) find all roots of a^3=d+lambda*b^3.
Skip a zero known systematic group or zero target in these
two cases.

For (c,d), skip c=0. Define k=lambda*c^3, which is nonzero.
The quantity A=a^3 satisfies
A^2+d*A+k=0, because A+lambda*b^3=d and
A*lambda*b^3=k. Find both possible A in K as follows.
At d=0 the square root of k is one possible value.
At d!=0 divide by d^2 to solve w^2+w=k/d^2; if a solution
exists, A=d*w and A+d are all solutions. For each A
that is a cube, choose a base cube root x and let y=c/x.
Both are nonzero because k!=0 and hence A!=0. Then
lambda*y^3=k/A=d+A, so (x,y) meets the stated check.
Rotate by (zeta*x,zeta^2*y) and (zeta^2*x,zeta*y) to get
all roots with that A and the same c. The backward
reasoning shows every good compatible message is generated.

The macro for cube roots uses e=(2h+1)/3, computed exactly.
For a cubic z!=0, z^h=1 and (z^e)^3=z. The same identity
holds at zero. A proposed root x=z^e is checked by cubing.
If valid, its multiples by zeta and zeta^2 give all
roots by degree, with duplicates only at zero. In the
charts admitting cube roots the target is nonzero.
For its fixed addition chain use
e=1+56*sum_(i=0..82) 2^(6i).
Compute z^7 as z^2*z^4*z, raise to 2^3 to get z^56,
then take its geometric power of length 83 with step 6
and multiply by z. For inversion of x!=0 use
x^-1=(x^(2^499-1))^2. This takes geometric length 499
with step one, then one squaring. To compute a geometric
power of length L of input T, set G=T and a=1.
For each remaining binary digit of L after its leading
one, replace G by G*G^(2^(step*a)). If the digit is one,
also multiply by T^(2^(2*step*a)). Update the compile-time
a to 2a plus the digit. The identities follow by splitting
the geometric exponent into adjacent intervals.

The linear map taking z to a possible solution of
w^2+w=z is AS_map(). The code computes a basis of the
image of x^2+x with tracked preimages, using descending
pivot positions. Reducing z against the fixed pivots
and accumulating their corresponding preimages is linear
over F2. On the image the residual is zero, and its
output is a preimage. The program verifies the equation
after this map, so for z outside the image it discards
its result. The two preimages for an attained z differ
by one, since the kernel is {0,1}. The declared computation
checks the AS map on the image basis, the exponent
identities and the chart rules on sample messages.

## Fixed RAM representation and macros

A field element occupies eight consecutive words, low first,
with 52 live bits in the top word. The memory arrays BASE
in construction.py start at 48 and are consecutive eight-word
arrays in ARRAY_LIST order. NEXT and LAST each contain
109 words starting at 312 and 421. MULT is the sixteen-row
scratch table at 530 with eight words per row. Fixed reduction
tables for shift amounts one and four occupy addresses 658
and 660, with two and sixteen entries respectively. An entry
at c is the carryless product of c with P_low=P+X^500, which
fits a word. For each of the 27 fixed linear maps in map_tags(),
starting at 676, store 63 tables of 256 entries, eight words
each. At byte position j and input byte d, the entry is the
map applied to the polynomial d*X^(8j), with its nonexistent
input bits above position 499 treated as zero. This is
enumerated by map_words(). The maps are the Frobenius powers
needed by the two chains and quadratic square root,
constant multiplications by lambda, its inverse, zeta and its
square, and the AS map. Addresses through the exclusive extent
3484324 fit the model, and the working arrays remain disjoint
from all fixed ranges.

All loops over literal offsets, table positions and chart
number in the macros below are code generation. All field
macro calls inline their primitive instructions. The only
runtime control is the stated chart guards and two-way
branches, plus the traversal using labels. Array arguments
of a macro stand for literal addresses at a particular
expansion. Moves and literal constants are free. Array
COPY takes two primitives per word, load and store.
Array XOR takes two loads, xor and store per word.
Thus COPY costs 16 and XOR costs 32. ZERO(A) ORs all
eight loaded words into an initially zero register and
compares with zero, costing 17. EQ(A,B) similarly loads
both words, XORs them and ORs into an accumulator, then
compares with zero, costing 33. A branch using one of
these flags costs one additional primitive.

SHIFT_s multiplies the field registers R by X^s modulo P
for s=1 or 4, emitting

    c := shr(R[7],52-s)
    for j=7 down to 1:
        h := shr(R[j-1],64-s)
        R[j] := shl(R[j],s)
        R[j] := or(R[j],h)
    R[0] := shl(R[0],s)
    R[7] := and(R[7],2^52-1)
    addr := add(REDSHIFT_s,c)
    h := load(addr)
    R[0] := xor(R[0],h)

The descending order retains the bits for carries.
Replacing the high coefficients by c*P_low performs
the modular reduction entirely in the low word.
This costs 1+21+1+1+3=27 operations.

APPLY(F,A,B) evaluates one of the fixed linear maps.
Initialize eight accumulators to zero. For each byte
j=0..62 of A use a literal load, shift and AND to get
the byte, multiply by eight and add the byte table base.
For each output word index add its literal value to
that address, load the word and XOR it into its
accumulator. Store the eight sums into B at the end.
Linearity gives the specified value. The final byte's
unused bits are zero. The count is 63*(5+24)+8=1835.

MUL(A,B,D) loads A into eight registers and stores
zero and A to rows 0 and 1 of MULT. Apply SHIFT_1 three
times, storing the intermediate vectors in rows 2,4,8.
For the eleven remaining nonzero rows, in increasing
order, use two loads and XOR then store each word of
the two earlier rows indexed by the least set bit and
the remaining bits. Thus row d is the product of A and
the degree-at-most-three polynomial d. Move R to zero.
For each literal nibble of B from high to low, load,
shift and mask it, multiply by eight, and add MULT.
Apply SHIFT_4 to R, then add each literal word index to
the row address, load and XOR. Horner's rule gives the
product. Store all eight words after the last nibble.
This uses
8+16+3*(27+8)+11*32+125*(5+27+24)+8=7489 operations.
Both APPLY and MUL finish all operand reads before
writing their destination so the stated in-place calls
are valid.

CUBE(A,B) is APPLY of the first Frobenius map into scratch
followed by MUL of that value with A into B. It costs 9324.
INV(A,B) copies A into G and uses the length-499
geometric chain in the preceding section. Each doubling
or extra term uses one APPLY and one MUL. There are eight
doublings and six extra terms. Use one last APPLY for
squaring. The cost is 16+14*7489+15*1835=132387.
ROOT(A,B) computes one candidate cube root in B and
returns its equality flag. Emit two APPLY to get z^2,z^4,
two MUL to form z^7, and one APPLY for raising to 2^3.
Copy the resulting T to G. The fixed length-83 chain
takes six doublings and three extra terms. Multiply G
by the original z to write B. Use CUBE on B and compare
with the original value by EQ. These are thirteen APPLY
and thirteen MUL including the cube check, one COPY and
one EQ, costing 121261. Calls to these macros use disjoint
scratch from the supplied input and output operands, with
G, H, J and the other temporary arrays in BASE reused after
each macro finishes.

## Chart evaluation in the word RAM

Store the original received field groups as Y0..Y3. In
each VISIT copy KEY0 and KEY1 to the corresponding candidate
groups V0..V3 for its chart, in 32 operations. During chart
evaluation the V groups describe the candidate field values
a,b,c,d, and used fields are restored to their required
values before every SCORE. SCORE computes distance by
loading each of the 32 V and Y word pairs, XORing,
popcounting and adding to Dist, costing 160. While
scanning the sixteen a,b words in low-to-high order,
also load the corresponding BEST word and apply ltu
and eq, then select the newly compared ordering when
not equal. This costs 64. Compare Dist with BestD by ltu
and eq, AND equality with the smaller-message flag, OR
with the strict-improvement flag and select BestD, at
five operations. For each systematic word load the
candidate and stored best word, select and store the
result at BEST0 or BEST1, at another 64 operations.
Thus SCORE costs 293 and implements precisely the numeric
tie rule. Field padding zeros make the order of the word
comparison agree with the order of a+2^500 b.

The chart bodies emit the following after the copies.
Each test for zero or equality is the named macro and a
conditional branch to EXTEND on an invalid key, or to the
next unrolled root on an invalid root. Except where noted,
completion falls through to EXTEND.

- For (a,b), guard both fields by ZERO. Compute V2 with
  MUL(a,b) and V3 as CUBE(a) XOR the lambda map of CUBE(b).
  Then SCORE. The guard costs 36, arithmetic 28004, SCORE
  293 and copies 32, below 30000.
- For (a,c) or (b,c), guard the two known fields, invert
  the known systematic field with INV and MUL by c to
  form the other. Emit two CUBE, APPLY(lambda), XOR to
  form d, then SCORE. The count is at most
  32+36+132387+7489+20515+293, below 165000.
- For (a,d) or (b,d), guard the known systematic field,
  form TARGET for the unknown cube using one CUBE, one
  APPLY of lambda or lambda_inv and one XOR, and guard
  TARGET. Emit ROOT(TARGET,ROOT) and test its flag.
  Evaluate ROOT, then its zeta and zeta2 maps, in
  three unrolled steps. The first copies the root to
  the missing candidate group, the others use APPLY to
  write it. In each step form c by MUL(a,b) and SCORE.
  Cost is at most
   32+36+9324+1835+32+121261+1+16+2*1835+3*(7489+293),
  below 165000.
- For (c,d), guard c and compute K=lambda*CUBE(c).
  Test d with ZERO and branch to the zero case. In the
  nonzero case emit INV(d,INV), APPLY(1,INV,INV2),
  MUL(K,INV2,TARGET), APPLY(AS,TARGET,ASROOT),
  APPLY(1,ASROOT,TEMP), XOR(TEMP,ASROOT,TEMP), and
  EQ(TEMP,TARGET), branching to EXTEND on failure.
  MUL(d,ASROOT,A0) and XOR(A0,d,A1), then jump over
  the zero case. In the zero case APPLY(499,K,A0)
  and COPY(A0,A1). This gives two root slots, with
  duplicate processing of the sole root allowed.
  For each slot emit ROOT to BASEA and branch to the
  next slot or EXTEND on failure. Emit INV(BASEA,INVA)
  and MUL(c,INVA,BASEB). COPY the two bases to V0,V1,
  then SCORE. Twice apply the opposite zeta maps to
  the bases to set V0,V1 and SCORE each rotated result.
  This emits at most two valid-root bodies. The
  prelude, including copies, guards, both alternatives'
  listed costs as an upper bound and their branches,
  is below 172000. Each root body is at most
  121261+1+132387+7489+32+4*1835+3*293,
  below 270000. The combined bound of 725000 covers
  the prelude and both bodies.

All field arrays in the chart proof are the fixed scratch
ranges from construction.py. The scratch results of two cubes
in an XOR are placed in C1 and C2 before combining, with the
lambda map done in place. Word-by-word XOR stores its result
after reading both operands at that index, so its uses in
place preserve correctness. In the inverse and root chains
the application scratch G,H,J is disjoint from the supplied
input and output base arrays. The field macros require at
most sixteen simultaneous field registers. Fewer than
twenty additional arithmetic temporaries suffice since only
the current word, address, and a fixed number of scalar
flags or shifts are live at each literal step. Depth and
BestD add two persistent registers. This is within 64.

## Traversal, output and total operations

The compiler extracts the four received groups into Y0..Y3.
A field piece of width at most 64 at fixed bit offset s
within input word j is obtained by loading that word and
shifting right by s if nonzero. On crossing a boundary,
load word j+1, shift it left by 64-s and OR it in.
Mask to the piece width if it is below 64. This takes at
most six primitives and extracts zero-padded pieces.
Each of 32 pieces is stored to its field array. This costs
at most 224 in total. Initialize BestD=2001 by a free move.

For each of six literal chart pairs emit two COPY from its
received groups into KEY0 and KEY1, set depth to zero and
store zero at NEXT[0]. This is 33 operations. Each chart
runs the following traversal, with labels local to the
chart and DONE falling into the next chart's initialization
or the final output.

    VISIT:
        emit this chart's VISIT body
    EXTEND:
        flag := eq(depth,108)
        branch if flag to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        flag := eq(i,1000)
        branch if flag to BACK
        j := add(i,1)
        store(adr,j)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        FLIP(i)
        jump to VISIT
    BACK:
        flag := eq(depth,0)
        branch if flag to DONE
        adr := add(LAST,depth)
        i := load(adr)
        depth := sub(depth,1)
        FLIP(i)
        jump to EXTEND

A branch not taken falls through. The mask enumeration visits
the increasing subsets of coordinates up to depth 108, once
each including the empty subset. Each node advances its own
NEXT before entering a child, which starts at the next
coordinate. LAST records the flip that is undone on returning
to the parent. FLIP uses four primitives to choose the
field base and its within-field index, then the seven
ordinary operations for a single bit update:

    first := ltu(i,500)
    j := sub(i,500)
    j := select(first,i,j)
    base := select(first,KEY0,KEY1)
    word := shr(j,6)
    off := and(j,63)
    mask := shl(1,off)
    addr := add(base,word)
    v := load(addr)
    v := xor(v,mask)
    store(addr,v)

The subtraction wraps as a word if i<500 but its result
is unused in that case. The selected index is always
0..499, so no padding is flipped. At an entering edge
the 21 ordinary operations with a seven-primitive flip
become 25. On each return the corresponding 13 become
17. A failed extension costs at most six and root
termination two. If N=sum_(i=0..108) binom(1000,i),
there are N visits, N-1 entries and nonroot returns,
N failures and one root termination, giving traversal
overhead at most 48N per chart.

At final output load the eight pieces of each best
systematic field, pack them into consecutive intervals
starting at bits zero and 500, and store to the required
output words. For each piece at zero word offset, store
it directly after loading. Otherwise shift it left,
load the output word, OR and store. On a crossing also
shift the original piece right by 64 minus the offset,
load the next word, OR and store. Each piece takes at
most nine operations. Intervals are disjoint and filled
in increasing order, so a direct store is only into a
word with no earlier assigned bits, and OR packing
preserves all earlier bits. Initially zero output words
and the padded field pieces give zero final padding.
The program halts after packing. Output cost is at most
144. The total fixed initialization and finish cost is
at most 224+6*33+144=566, hence is bounded by 600.

Summing the six VISIT bounds and all traversals gives
the all-input bound W=1415288*N+600. The declared
experiment checks W^1000<=2^509851 using exact
integers. The complete chart decoder and operation
count use the listed primitive program and deterministic
tables. The declared experiments verify the numerical
certificates and sample chart identities.
