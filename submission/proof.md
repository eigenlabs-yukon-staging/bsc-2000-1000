# Three incremental pair searches for a two-coordinate field code

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
encoded groups within that pair's total distance cap. The
caps in lexicographic pair order (a,b), (a,c), (a,d), (b,c),
(b,d), (c,d) are (111,111,109,111,109,108). Write R_I=111
for the systematic pair, R_C=111 for mixed pairs involving c,
R_D=109 for mixed pairs involving d, and R_P=108 for the
parity pair. Among all eligible messages
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
of weight at most t with some pair of groups with their combined
weight within that pair's cap. Let V_t=|E_t|. For a received vector y let
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
If u+v<=R_I, the permitted parity error weights p_0,p_1
satisfy simply p_0+p_1<=ell. Otherwise the condition is

    p_0+p_1<=ell and
    (p_0+p_1<=R_P or p_0<=R_C-min(u,v)
                  or p_1<=R_D-min(u,v)).

For the latter, the cheapest mixed partner among the
systematic groups has weight min(u,v), separately for either
parity group. The parity sets are invariant under permutations
within each 500-bit group. Their volume J sums
binom(m,p_0)binom(m,p_1) over the specified weight pairs.

Let D=42. A half-mask of weight between D and m-D inclusive
is balanced. For that half, its signed Hamming-shell
coefficient at weight w has absolute value at most
C_w=[x^w](1+x^2)^D(1+x)^(m-2D). In fact, its exact
signed sum is the coefficient of (1-x)^s(1+x)^(m-s) for
mask weight s. Pair at least D opposite factors and take
absolute values coefficientwise. Surplus pairs are bounded
by (1+x)^2 in place of 1+x^2. For the other half the
trivial bound is its shell cardinality. If the first half is balanced,
the absolute signed sum over the parity set is at most the
sum of C_(p_0)*binom(m,p_1) over the permitted pairs. If
the second half is balanced, use the corresponding sum
of binom(m,p_0)*C_(p_1). Take C to be the maximum of
these two computed bounds. Translating
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
ball up to ell or the one determined by the ordered pair
(R_C-min(u,v), R_D-min(u,v)), with each threshold
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
At parity total w<=R_P all the shell is included. When
w>R_P let the clamped thresholds be A and B. The shell
includes p_0<=A or p_0>=w-B. Its prefix/suffix sum counts
the union exactly, taking the whole shell if the intervals
meet or are adjacent. This works for J and for each
orientation of the coefficient bound.

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
13793970757048787/250000000000000000. Normalization is
by all Q messages, including the ones the decoder excludes.

## The six field charts

The chart for a specified pair of group values produces
exactly its compatible (a,b) with a,b nonzero. Chart pairs
are visited in lexicographic pair order. On each chart
enumerate all two-group received perturbations of weight
at most that pair's cap. Thus the union of the output messages is exactly
the eligible set in the decoder definition. Duplicate messages
do not alter the minimum or numerical tie rule.

Pairs (a,b), (a,c) and (b,c) determine the message directly:
skip keys with a required nonzero value zero, and for
(a,c) or (b,c) recover the other component by division.
For (a,d) find all roots of b^3=lambda^-1*(d+a^3).
For (b,d) find all roots of a^3=d+lambda*b^3.
Skip a zero known systematic group or zero target in these
two cases.

For (c,d), skip c=0 or d=0. A good message has both
c!=0 and d!=0, because a^3+lambda*b^3=0 would make lambda
a cube. Define k=lambda*c^3. The two cubic quantities
a^3 and lambda*b^3 must be the roots of
A^2+d*A+k=0, since their sum is d and their product is k.
Compute a root w of w^2+w=k/d^2, or return no candidate
if none exists. Put A=d*w. It is nonzero and the other
root is A+d. One need only raise this one A to the cube
root exponent e below. Write t=A^e. If t^3=A, use
a=t, b=c/a. If t^3=zeta^2*A, then
(A^h)^2=zeta^2, so A^h=zeta. In this case let
b=t*lambda^(-e), whose cube is
A*zeta^2/(lambda*zeta^2)=A/lambda, and set a=c/b.
In the remaining case there is no good solution. Indeed,
for any good solution A must be either the cube a^3 or
lambda times the cube b^3, and those two cases force the
two stated checks on (A^e)^3. Conversely, in the first
case lambda*b^3=k/A=A+d, and in the second
a^3=k/A=A+d, so the equations hold. The cases are
disjoint since lambda is a noncube. Rotate the found base
(a,b) by (zeta*a,zeta^2*b) and (zeta^2*a,zeta*b).
These three pairs exhaust all roots of the chosen cubic
and their corresponding division values.

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
arrays in ARRAY_LIST order. NEXT, LAST, INEXT and ILAST each
contain 112 words and start at 368,480,592,704 respectively.
MULT is the sixteen-row scratch table at 816 with eight words
per row. Fixed reduction tables for shift amounts one and
four occupy addresses 944 and 946, with two and sixteen entries.
An entry at c is the carryless product of c with P_low=P+X^500.
For each of the 27 fixed linear maps in map_tags(), starting
at 962, store 63 tables of 256 entries, eight words each.
At byte position j and input byte d, the entry is the map
applied to the polynomial d*X^(8j), with input bits above
position 499 treated as zero. This is enumerated by map_words().
The maps are the Frobenius powers needed by the two chains,
constant multiplication by lambda, lambda^-1, lambda^(-e),
zeta and its square, and the AS map. Starting at address 3484610, DIRECTIONS
has 500 rows of eight scratch words. SAVED starts at
3488610 and has 112 rows of sixteen scratch words. AFFINE
starts at 3490402 and has room for 500 direction tables in scratch of
63 byte slices, with the same byte layout as APPLY.
The exclusive memory extent is 68002402, and all
ranges are disjoint.

All loops over literal offsets, table positions and chart
number in the macros below are code generation. All field
macro calls inline their primitive instructions. Runtime control consists of chart guards and branches,
the affine-table generation loops, and the traversal. Array arguments
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
ROOT_POWER(A,B) computes the candidate root A^e in B
and its cube in CHECK. ROOT(A,B) uses this and then returns
the equality flag with the original value A. Emit two APPLY to get z^2,z^4,
two MUL to form z^7, and one APPLY for raising to 2^3.
Copy the resulting T to G. The fixed length-83 chain
takes six doublings and three extra terms. Multiply G
by the original z to write B. Use CUBE on B into CHECK.
ROOT_POWER uses thirteen APPLY and thirteen MUL, including
the final cube, and one COPY, costing 121228. ROOT adds
one EQ for a total of 121261. Calls to these macros use disjoint
scratch from the supplied input and output operands, with
G, H, J and the other temporary arrays in BASE reused after
each macro finishes.

## Outer and inner chart arithmetic

Store the original received field groups as Y0..Y3. A chart
splits its candidate masks by group. It enumerates an outer
500-bit field and, for each outer key, all eligible inner
500-bit fields. The outer group is the first named in the
pair except for (c,d), where it is d. The total of outer and
inner mask weights is bounded by that chart's cap. For every
outer visit, copy KEY0 to its corresponding field among
V0..V3. The inner visit similarly copies KEY1 to its own V
field. Other V arrays are used for the recovered candidate
values and are set as needed before every SCORE. SCORE computes distance by
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

The outer visit begins with the COPY described above and a
ZERO guard on its known field, branching to the outer extension
on zero. Every good message has that field nonzero. Emit the
following precomputation once at this outer key.

- In (a,b), (a,c) and (b,c), build the direction tables
  and initialize the recovered fields as specified below.
- In (a,d), write a^3 to KNOWN using CUBE.
- In (b,d), write lambda*b^3 to KNOWN.
- In (c,d), write d^-1 to INV and square it into INV2.

The three charts without direction tables have outer visit
cost below 200000, including the initial COPY and guard.
The uniform bound over all six charts is 500000000, including
the generation described below. The immutable outer field and precomputed
arrays INV, INV2 and KNOWN are kept throughout the inner search.
Each inner evaluation begins with the 16-operation COPY of KEY1
and uses these bodies. A failed guard branches to the inner EXTEND.

- For (a,b), guard b and SCORE. For (a,c) or (b,c), guard c
  and SCORE. In all three charts the traversal maintains the
  recovered groups at their proper values including when the
  inner group is zero. A nonzero outer a or b means a nonzero c
  determines a nonzero recovered systematic group.
- In (a,d) form TARGET by XOR(d,KNOWN) and APPLY of
  lambda^-1. For (b,d) the target is XOR(d,KNOWN).
  Guard zero TARGET. Use ROOT(TARGET,ROOT), discarding
  a failed equality. For each of the base root and its two
  zeta rotations, write the missing systematic field, form
  c by MUL(a,b) and SCORE. The first uses COPY and the next
  two use APPLY. A bound covering either chart is
  16+32+1835+18+121261+1+16+2*1835+3*(7489+293),
  below 152000.
- In (c,d), guard c and form K=lambda*CUBE(c).
  Multiply K by precomputed INV2 to get TARGET.
  APPLY the AS map to TARGET into ASROOT, square ASROOT
  into TEMP and XOR ASROOT into TEMP. Use EQ to check that
  it matches TARGET and branch to EXTEND on failure.
  MUL(d,ASROOT,A0) forms A. Emit ROOT_POWER(A0,ROOT)
  whose CHECK is the cube of the proposed field root.
  Use EQ(CHECK,A0) and branch to the cubic orientation
  if true. Otherwise use APPLY of zeta2 to A0 into TEMP
  and EQ(CHECK,TEMP). Branch to EXTEND on failure.
  This orientation sets BASEB by APPLY of lambda^(-e)
  to ROOT, then INV(BASEB,INVA) and MUL(c,INVA,BASEA).
  Jump to the shared scoring code. In the cubic orientation
  copy ROOT to BASEA and emit INV(BASEA,INVA) and
  MUL(c,INVA,BASEB). It falls through to the shared code.
  SCORE the base pair with a COPY of each base into V0,V1.
  For each of the two rotations, apply two zeta maps from
  the base arrays to V0,V1 and SCORE. This covers all
  field chart candidates. For a successful noncubic
  orientation the count is at most
  16+18+11159+7489+1835+1835+32+34+7489+121228+
  34+1835+34+1835+132387+7489+1+32+4*1835+3*293
  =303001. The cubic orientation does not execute
  the extra equality, maps and jump. Failures skip
  some of these operations. The bound 310000 covers
  all inner visits in (c,d).

All field arrays are the fixed scratch from construction.py.
CUBE results are placed in C1 or C2 before the following
XOR or constant map. Wordwise XOR writes each word after
reading both of its operands at that index, allowing its
stated in-place use. Field macro scratch G,H,J and the other
arrays stated above are disjoint from the named input and
output operands at their calls. The field macros require at
most sixteen simultaneous field registers. Fewer than twenty
additional arithmetic temporaries suffice because only the current
word, address, and fixed scalar shifts or flags live at each
literal step. The two depths, the inner cap, and BestD add
four persistent registers, within 64.

### Incremental inner traversals for (a,b), (a,c) and (b,c)

Write u for the inner field and Z for the recovered field
other than d. In (a,b), u=b and Z=c=a*u. The parameters
S=a, alpha=lambda and D0=a^3 yield Z=S*u and
d=D0+alpha*u^3. In (a,c), u=c and Z=b; take S=a^-1,
alpha=lambda*S^3 and D0=a^3. In (b,c), u=c and Z=a;
take S=b^-1, alpha=S^3 and D0=lambda*b^3.
At the outer visit, compute S, alpha and D0 by the stated
COPY, INV, CUBE and APPLY macros, storing them in INV,
ALPHA and KNOWN. For (a,b) use a COPY for S and eight
literal stores for alpha. Initialize the Z array with
MUL(S,Yinner) and initialize V3 with CUBE(Yinner),
MUL(ALPHA, cube), and XOR with D0. The other known fields
are populated by the usual outer and inner visits.
This initialization uses at most one INV, three CUBE,
two APPLY, two MUL, four COPY, one XOR, eight literal
stores and the outer guard. It costs less than 200000.
Load S into eight registers. For each literal i from 0
through 499, store the registers at row i of DIRECTIONS
and apply SHIFT_1 except at the last row. Thus row i is
S*X^i. This costs 8+500*8+499*27=17481.

For each outer visit in these charts, regenerate AFFINE
using its alpha. For direction e=X^i, the table represents
the affine map alpha*(e*u^2+e^2*u+e^3). Its ith direction
is at AFFINE+i*AFFINE_STRIDE. Here is a counted
construction of the complete table.

1. Initialize F1,F2,F3 with copies of ALPHA. Maintain
   F1=alpha*X^(2*i), F2=alpha*X^i, F3=alpha*X^(3*i).
   For each direction i, copy F1,F2 to LIN,QUAD and
   start the byte loop at j=0. Use a cursor initialized
   at AFFINE before the direction loop and advance it
   by 2048 after each byte.
2. In byte j store entry zero as F3 when j=0 and as
   eight zero words otherwise. For each literal bit
   r=0..7 use column LIN XOR QUAD. For r>=4 with j=62,
   use the zero column instead. Obtain each of its eight
   words by two loads and XOR, with one optional select
   per word using the comparison j==62. In increasing
   order of the literal bit r, for each literal k in
   0..2^r-1 and each literal word index, add the
   entry and word offset to the cursor, load the
   old word, XOR the indicated column word, add the
   new offset and store at entry k+2^r. Each word
   can keep its column in a register while running the
   literal k steps. Then multiply LIN by X and QUAD
   by X^2 using SHIFT_1 respectively once and twice,
   with eight loads and eight stores to each array.
   These shifts advance the actual polynomial columns
   even for the four padding positions.
3. After its 63 bytes advance F1,F2,F3 by SHIFT_1
   respectively twice, once and three times. Store the
   updated arrays and advance i. Use increasing counters
   with literal end values 500 and 63. Emit comparison,
   branch, increment and jump for each loop, and pass
   through the byte action only for j=0 through 62.

The invariant at the column for index l=8*j+r is
LIN=alpha*e^2*X^l and QUAD=alpha*e*X^(2*l). The XOR and
entry expansion therefore produce the table specified
by derivative_words(i,alpha), including the constant
exactly once in each APPLY. Initializing row zero uses
at most 40 primitives including the test and either
literal zero stores or a copy with dynamic stores.
For the columns and entries each byte costs at most
8*(8*4+3*27+32) + 255*8*5. Even adding 100 for
row zero and byte loop overhead, this is below 14000.
Each direction's extra work, including two copies, six
shifts, loads and stores of the three arrays, and
direction loop overhead, is below 300. Thus table
generation costs less than
500*(63*14000+300)+100, below 442000000,
and the other precomputations fit within the full
500000000 cap. In building the table the outer depth
and BestD are retained, and at most eight shift registers
and twenty additional registers suffice. The table and
all named scratch ranges are disjoint.

Before the inner FLIP on descent save the old Z and V3
in the row SAVED+16*depth of the new depth. Emit one mul,
one add, and for each of sixteen literal offsets emit
add, load and store, costing 50. Form the address
DIRECTIONS+8*i using one mul and one add and XOR its row
into Z, emitting an add, two loads, XOR and store at
each of eight literal offsets. This costs 42 and changes
Z from S*u to S*(u+e).

For V3 use the affine difference alpha*(e*u^2+e^2*u+e^3).
Emit two operations to get base AFFINE+i*AFFINE_STRIDE,
then use APPLY with eight accumulators initially loaded
from V3 instead of zero. At each byte build its address by
adding the literal byte offset to the register base, then
adding eight times the byte. This adds one operation per
byte to APPLY. The total is 2+8+1835+63=1908.
Its input is the old KEY1. At the end V3 is the parity
value for u+e. Execute the ordinary FLIP immediately
after this update before entering the child visit.

The inner BACK in these charts restores Z and V3 from the
current depth row before decrementing depth. Emit mul and
add to form the row base, then per literal offset add the
offset, load the saved word and store it to the appropriate
array. The 50-operation restore and the usual FLIP recover
the parent. The zero guard and SCORE do not modify the
recovered arrays. This extension and restoration add
50+42+1908+50=2050 per nonroot inner node, pairing
its entry and exit. A visit itself costs 16+18+293=327,
so this fits in 2400 operations per inner node including
the root. The standard traversal work is accounted for
separately. The extra address and word temporaries reuse
the APPLY and word-copy registers, within the stated limit.

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

For each of six literal chart pairs, copy its received outer
field into KEY0, set outer depth to zero and store zero at
NEXT[0], costing 17. The outer search traverses 500 positions
up to the chart cap. Its VISIT executes the chart precomputation.
If a guard fails, it branches to its extension. Otherwise copy
the inner received field into KEY1, set inner depth zero,
compute the inner cap as chart cap minus outer depth, and
store zero in INEXT[0]. This costs 18 operations.
Immediately follow with the compiled inner search, which uses
the chart's inner VISIT body and exits to outer EXTEND.
The outer search then advances normally.

For either traversal the compiler emits the following template,
with array addresses and labels substituted. KEY is KEY0
or KEY1. NEXT and LAST are their outer arrays or INEXT
and ILAST for the inner search. The two depth registers
are distinct. The outer cap is a literal and the inner
cap is the computed register value. DONE is outer EXTEND
for the inner copy, and the next chart initialization
or the final output for the outer copy.

    VISIT:
        execute the corresponding inlined VISIT action
    EXTEND:
        flag := eq(depth,cap)
        branch if flag to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        flag := eq(i,500)
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

A branch not taken falls through. The three incremental inner templates additionally emit the
descent maintenance and ascent restoration above. The outer VISIT contains the
inner template before its own EXTEND label, so inner root exit
uses the outer EXTEND and no call or return is implicit.
The mask traversal advances the parent NEXT before entering a
child, starts the child at the next greater coordinate, and
records the flip in LAST for undo. Thus a cap-r traversal
visits exactly the subsets of size at most r, including the
empty subset, unless an outer guard skips the whole inner
search. FLIP emits seven primitives:

    word := shr(i,6)
    off := and(i,63)
    mask := shl(1,off)
    addr := add(KEY,word)
    v := load(addr)
    v := xor(v,mask)
    store(addr,v)

The index is always in 0..499, preserving padding. A successful
extension costs 21, a nonroot BACK costs 13, and a failed
extension costs at most six. Root BACK costs two, so the
traversal overhead is at most 40 per node summed over either
complete traversal. The outer visited node count is
O_r=sum_(i=0..r) binom(500,i). The total potential inner
node count across the outer tree is exactly
N_r=sum_(i=0..r) binom(1000,i), by choosing the two
disjoint parts of a subset with total size at most r.
Discarding an outer key can only decrease the count.

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
at most 224+6*17+144=470, hence is bounded by 600.

For a chart with radius r and inner evaluation bound C,
combine the costs as (C+40)*N_r+500000058*O_r. The outer
constant includes its 500000000 precomputation cap, 18 to
start an inner search, and 40 of traversal overhead.
Sum these six expressions with the chart costs and
radii in construction.py and add 600 for fixed work.
This is the stated integer W for every input and
every control path. The declared experiment checks
W^125<=2^64012 with exact integers. Correctness and
the worst-case bound follow from the listed primitive
program and fixed tables. The exact computation also
checks example field chart identities without executing
the full candidate enumeration.
