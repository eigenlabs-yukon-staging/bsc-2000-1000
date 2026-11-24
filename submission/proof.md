# Fourteen charts for the quadratic-tower field code

## Code and candidate rule

Let F be the binary field with modulus
P=X^250+X^59+X^40+X^32+1. The checker verifies
X^(2^250)=X modulo P and gcd one with X^(2^125)+X
and X^(2^50)+X. An irreducible factor of degree j
divides X^(2^l)+X exactly if j divides l: Frobenius
has order j on its quotient by the multiplicative
group order and the root bound, and fixing X fixes
the whole quotient. Thus the checks exclude all
proper factor degrees. Set xi=X^191. The checked
binary trace of xi in F is one, whereas x^2+x
has trace zero by telescoping. Hence
K=F[theta]/(theta^2+theta+xi) is a field of size
q=2^500. Use the two consecutive 250-bit
coefficient blocks, low coefficient first, for
each K element. Let lambda=X in F, h=(q-1)/3,
zeta=lambda^h and e=(2*h+1)/3. The checker finds
zeta!=1 and zeta^3=1, so lambda is not a cube.

Encode a,b as (a,b,c,d)=(a,b,ab,a^3+lambda*b^3),
or 2000 bits in four groups. The systematic layout
gives all 2^1000 distinct messages. A field is good
when its two F coefficients are nonzero and distinct.
Only messages with both a and b good may be decoded.
They qualify by one of these chart caps.

- The combined noise weight of any of the six
  full group pairs is at most 105.
- For either full systematic group, combine its
  weight with one half of the other systematic
  group and one half of c. Any such subtotal
  at most 106 qualifies. There are eight choices.

Choose the qualifying message of minimum full
distance, breaking ties by numeric a+q*b. The set
is nonempty. A received group with both halves zero
can be changed to halves 1 and X in two flips. If
just one half is zero, change it to 1 or X distinct
from the other. If both halves are equal and
nonzero, a single flip by 1 or X in one of them
makes them distinct and nonzero. At most four flips
in a,b qualify via their ordinary pair.

## Parity differential and collisions

For a nonzero systematic difference (s,t), the
parity differential is affine over F2 with values

    at+bs+st,
    s*a^2+s^2*a+s^3+lambda*(t*b^2+t^2*b+t^3).

For s!=0, fixing the first value makes b affine
in a with slope t/s. The linear part of the second
is (s^3+lambda*t^3)/s^2 times a^2+s*a.
The coefficient is nonzero since lambda is not a
cube. The kernel is {0,s}. For s=0,t!=0 the first
value fixes a and the remaining b map has kernel
{0,t}. The joint parity difference is therefore
uniform on an affine hyperplane, with mass 2/Q
per point for Q=q^2 uniform messages.

The trace pairing on K is nonsingular: the
binary trace polynomial is nonzero by its degree
less than q, and multiplication by a nonzero value
is a bijection. Hence each binary parity mask is
uniquely Tr(alpha*c+beta*d). Its linear differential
annihilator equations are

    alpha*t + beta*s^2 + sqrt(beta*s) = 0,
    alpha*s + beta*lambda*t^2 + sqrt(beta*lambda*t) = 0.

They use trace invariance under squaring. When
beta=0, a nonzero mask admits only zero direction.
When alpha=0 and beta!=0, squaring gives two
nonzero quartics and at most sixteen solutions.
Otherwise square the first equation to obtain
t^2=A*s^4+B*s with A!=0. Its square substitutes
for t^4 in the squared second equation, giving
t=C*s^8+D*s^2 with C!=0. Compare its square to
the expression for t^2 to get a nonzero degree
sixteen equation for s. Then t is determined.
All cases include zero, so at most fifteen
nonzero directions share one mask.

Fix actual noise z and competing systematic noise
of weight j. XOR fixes the direction. Allow the
competitor any parity noise in the ordinary ball
of radius t-j, ignoring chart admission. Let J
be the ball count. If at least one 500-bit mask
half has at least D0=42 ones and at least D0 zeros, the signed ball
count is at most H in absolute value, with H the
sum through t-j of coefficients of
(1+X^2)^D0*(1+X)^(1000-2*D0).
Indeed each shell's signed sum is the coefficient
of (1-X)^w*(1+X)^(1000-w) for mask weight w.
Pair D0 factors of opposite signs and bound
by the triangle inequality. Then H<=J.
Translation by actual noise only changes the
sign. The hyperplane intersection has size at
most (J+H)/2. The number of directions lacking
such a balanced half is at most
Kbad=15*(2*sum_(i<D0) binom(500,i))^2.
Use the trivial intersection J on these.
For N=binom(1000,j) patterns, summing
N*(J+H)+min(N,Kbad)*(J-H) over j=0..t gives A_t.
Thus the expected competitor count for any
fixed z is bounded by A_t/Q.

## Exact admitted volume and success

For total noise through t<=240 let V_t count
errors meeting some chart. Write A,B,C,D for
full weights and v_X for the smaller half
weight of group X. Put s=A+B, p=C+D, m=min(A,B),
u=min(A+v_B,B+v_A). All charts fail exactly if
s>105, p>105, m+C>105, m+D>105 and u+v_C>106.
Thus only 106<=s<=134 and 106<=p<=240-s are
needed for rejected shells. For a full weight w,
define T(w,k)=sum_(j=k..w-k) binom(250,j)binom(250,w-j),
the count with both halves at least k, with the
whole count for k<=0 and zero above w/2.
At fixed A<=B the count for u>=l is
T(A,l-B)*T(B,l-A). Differences give equality
counts, doubled when A<B. The possible l range
is A..min(A+floor(B/2),B+floor(A/2)).
For fixed p,m,u the remaining failure count is
sum_C T(C,107-u)*binom(500,p-C) over C>=0,
p-C>=0, C>105-m and p-C>105-m.
Subtract rejected shell counts from binom(2000,t)
and cumulate to get V_t. The checker uses exactly
these integer sums.

There are G=((2^250-1)*(2^250-2))^2 good
messages. For a received y let m_y count admitted
good messages of distance at most t. Then
sum_y m_y=G*V_t. For each admitted true error,
the competitor estimate gives ordered distinct
collisions at most A_t after summing over all Q
messages, so sum_y m_y^2<=V_t*(G+A_t).
Let F_t count received y with m_y>0, divided by Q.
Cauchy-Schwarz gives
F_t>=V_t*G^2/(Q*(G+A_t)). Also for integer l>=1,
(x-l)(x-l-1)>=0 for positive integer x implies
l(l+1)>=(2*l+1)x-x^2, valid with zero left side
when x=0. Summing yields
F_t>=V_t*(2*l*G-A_t)/(Q*l*(l+1)).
Use l=1+floor(A_t/G), floor the two lower bounds
and maximize with the previous D_(t-1), starting
at D_(-1)=0. Thus D_t<=F_t and D_t is nondecreasing.

At winning distance t, each received y contributes
Q^-1*9^(2000-t)/10^2000 to average decoding
success. The winning-distance cumulative count
divided by Q is F_t. Hence
sum_(t<=240)(D_t-D_(t-1))*9^(2000-t)/10^2000
lower-bounds success, by summation by parts
with nonnegative coefficients. The declared
computation obtains error at most
11166727434631111/200000000000000000.

## Field chart solutions

In (a,b) both values are known. In (a,c) and
(b,c), divide c by the known nonzero factor.
For (a,d) solve b^3=lambda^-1*(d+a^3).
For (b,d) solve a^3=d+lambda*b^3. For a nonzero
target z, test z^e by cubing. If it is a cube,
z^h=1 and (z^e)^3=z, and all roots are the tested
root and its two zeta rotations.

In (c,d), discard zero c or d. Good messages
never have either zero since lambda is not a
cube. Set k=lambda*c^3. The two quantities
a^3 and lambda*b^3 are roots of A^2+d*A+k=0.
Use the linear AS map for a possible root of
w^2+w=k/d^2 and test the equation. AS_map
eliminates against a fixed basis of 499 image
vectors of x^2+x, with their preimages tracked.
On the image it returns a preimage, and outside
it the test rejects. The roots differ by one.
Put A=d*w, which is nonzero. Test R=A^e. If
R^3=A, take a=R and b=c/a. If R^3=zeta^2*A,
then A^h=zeta from 3e=2h+1, and b=R*lambda^-e
has cube A/lambda, so take a=c/b. Otherwise
neither required cubic orientation holds.
In either passing case the other cubic value
is k/A=A+d. The two base-pair rotations by
(zeta,zeta^2) and (zeta^2,zeta) exhaust the
solutions. Retain only good a,b.

To compute z^e use
e=1+56*sum_(i=0..82) 2^(6*i).
Form z^7 from z^2,z^4,z, then z^56 by
Frobenius step three and the length-83
geometric power with step six, multiplying
by z at the end. For x^-1 square the length-499
geometric power with step one. For a geometric
power of T, start G=T and a=1. As each binary
digit of the length follows its leading one,
replace G by G*G^(2^(step*a)); if the digit is
one also multiply by T^(2^(2*step*a)). Update a
to twice a plus that digit. Splitting the exponent
sum into adjacent ranges proves this chain.

For a split chart v is the known full systematic
value and u the other. Multiplication by v has
matrix [[v0,xi*v1],[v1,v0+v1]] on F coefficients.
All entries are nonzero for good v. If p0 is a
known coefficient of u and q0 a known coefficient
of c, use their row's constants t0,t1 for the
known and unknown columns. Recover the unknown
coefficient as t1^-1*(q0+t0*p0). The fourth
group d is Dbase+alpha*u^3, where Dbase=v^3 and
alpha=lambda for v=a, or Dbase=lambda*v^3 and
alpha=1 for v=b. Retain good u only. Each chart
recovers exactly the codewords matching its
known coordinates, so enumerating the specified
known-bit radius and scoring all recovered
words implements the decoder.

## RAM representation and arithmetic

Store a field element as two consecutive four-word F
halves. The order is low coefficients first, with six
zero high bits in the last word of each half. This
represents the same value as the two external bit groups
without the internal padding. An F element takes four
words and a K element takes eight. The 47 arrays in
ARRAY_LIST begin at address 48, with the BASE addresses
and eight words per array. The four traversal stacks
NEXT,LAST,INEXT,ILAST have 107 words each, beginning
at 424,531,638,745. The 16-row, four-word base product
scratch MULT begins at 852. The two and sixteen
reduction entries begin at 916 and 918. XI_MAP begins
at 934, with 32 slices of 256 four-word entries given
by xi_entry. MAP_BASE begins at 33702 and contains
the 29 maps of map_tags, each with 64 byte slices and
256 eight-word entries per slice. For both kinds of
map the table entry for a byte is the map value on
that byte in its named position, with all padding
coordinates contributing zero. The XI_MAP applies to
a single four-word base field. Each 64-slice map
uses 32 bytes in each half. The declared map_columns
and map_words functions completely determine this
fixed read-only range.

DIRECTIONS starts at 3834790 and has 500 rows of
sixteen words each. SAVED starts at 3842790 and has
107 rows of 24 words each. AFFINE starts at 3845358,
with 500 direction slots each having the same size
as one fixed 64-slice map. These are writable scratch,
reused across charts and keys. The exclusive end is
69381358, below 2^32, and the above ranges are disjoint.

Fixed index loops, literal cases and macro uses below
generate inlined instructions. The counted runtime
control consists of chart guards, traversal, base
inversion in the split setup if emitted with loop
control as described, and the affine table builder.
All field storage uses padded words, not operations on
wide integers. COPY and XOR use two and four
primitives per word, respectively, so have costs 16
and 32 for K and 8 and 16 for F. XOR can operate in
place, since it reads both words before writing each
index. ZERO on a K array uses eight loads, eight ORs
into an initially zero register and an equality test.
EQ on two K arrays similarly loads both words, XORs
each and ORs into the register, then tests it, for 33.
Add one branch operation to either test if branching.

For a word-aligned four-word base field, BSHIFT_s for
s=1 or 4 multiplies the registers R by X^s modulo P.
Emit

    c := shr(R[3],58-s)
    for j=3 down to 1:
        H := shr(R[j-1],64-s)
        R[j] := shl(R[j],s)
        R[j] := or(R[j],H)
    R[0] := shl(R[0],s)
    R[3] := and(R[3],2^58-1)
    address := add(REDSHIFT_s,c)
    H := load(address)
    R[0] := xor(R[0],H)

Each REDSHIFT entry is the carryless product of the
discarded c with P+X^250 by reduction_table.
Its degree is at most 62. The descending shifts
preserve the old words for carry extraction. The macro
costs 1+9+1+1+3=15.

BMUL(A,B,D) is the following fixed Horner computation.
Load the four words of A. Store zero and A to rows 0
and 1 of MULT. Build rows 2,4,8 by applying BSHIFT_1
to the register value three times, storing after
each shift. For the remaining eleven nonzero rows
in numeric order, XOR the earlier row for the least
set bit with the row for the remaining bits, using
two loads, xor and store per word. These rows are the
products of A with their nibble polynomials.
Move the four R registers to zero. For each of the
63 literal nibbles of B, from high to low, emit load,
shift, mask, multiply by four and add MULT to obtain
the row address. Apply BSHIFT_4 to R. For each word,
add the index to the address, load the row word and
XOR it into R. Store the four registers to D when
done. This is Horner's rule for multiplication, and
its count is
4+8+3*(15+4)+11*16+63*(5+15+12)+4=2265.
The input words are read before the destination is
written. MULT is separate from all operands.

An APPLY map on K uses eight zeroed accumulators.
Each byte uses literal load, shift and AND of the
appropriate half's word, then multiply by eight and
add its literal slice base. Per output word add its
index to this address, load and XOR into the
accumulator. Store the eight sums at the end.
This costs 64*(5+24)+8=1864. Base APPLY for XI_MAP
uses the same rule with four accumulators and 32
bytes, costing 32*(5+12)+4=548. The identity follows
by F2 linearity. The last byte of each F half has only
two live positions.

For MUL on K form A0+A1 and B0+B1 into the first
halves of M0 and M1. Use BMUL to put A0*B0 into
M2, A1*B1 into M3 and (A0+A1)*(B0+B1) into M4,
using the first half of each listed scratch array.
Apply XI_MAP to M3 in place. The two output halves
are M2+M3 and M4+M2, by the basis multiplication.
This costs 4*16+3*2265+548=7407. No chart operand
is in M0 through M4. The final writes can alias an
input because all input words were consumed earlier.

CUBE uses APPLY of squaring into X2 and then MUL with
the original operand, at 9271 operations. An invocation
whose input occupies X2 instead uses another disjoint
scratch array such as X4, because the original must
be available to MUL. The two geometric chains of the
ordinary charts use APPLY for each fixed Frobenius
power and MUL for each product. INV takes one COPY,
14 MUL and 15 APPLY, totaling 131674. ROOT_POWER
takes one COPY, 13 MUL and 13 APPLY including the
final cube written to CHECK, totaling 120539.
ROOT adds one EQ, giving 120572. In the geometric
chains use G for the running product and H,J for
term images, leaving the input and any other retained
chart value disjoint. ROOT_POWER uses X2, X4, T for
the first power before its chain, and can reuse X2
for the last cube once the root has been written to ROOT.
As with other field operations, each macro expands to
fixed reads and writes and may reuse inactive scratch.

A good-field guard loads each pair of corresponding
base words, XORs them, and ORs each word and their XOR
into three initially zero registers. This is six
operations per word. Use three ltu(0,accumulator) comparisons,
two ANDs of the resulting nonzero tests, and one branch.
The resulting count is 30. Use separate tests of a and b
before each SCORE. The score itself computes full
distance by loading each of 32 pairs of Y and V
words, XORing, popcounting, and adding to Dist at
160 primitives. During those distance loads, for each of the sixteen
a,b words in low-to-high order, retain that word in
its current register. Load the corresponding BEST
word, compare by ltu and eq, and update the smaller flag
using select whenever words differ. This costs 64
and implements the numerical tie rule because the word
order respects the coefficient order. Use ltu, eq, AND,
OR and select for the comparison and update of BestD.
For each systematic word load candidate and prior best,
select and store to BEST. This costs 64 more, giving
293. Good tests and SCORE together cost at most 353.

## Counted ordinary chart bodies

For each full pair the outer key is the first group,
except for (c,d), whose outer key is d. At an outer
visit copy the key to its V field and branch past a
zero value. Prepare the indicated arrays once.

| Pair | Retained outer quantities |
| --- | --- |
| (a,b), (a,d) | KNOWN=a^3 |
| (a,c) | KNOWN=a^3, INV=a^-1 |
| (b,c) | KNOWN=lambda*b^3, INV=b^-1 |
| (b,d) | KNOWN=lambda*b^3 |
| (c,d) | INV=d^-1, INV2=d^-2 |

One COPY, guard, INV, CUBE and constant APPLY cost
below 200000, covering every row. The inner visit
copies KEY1 to its matching V field at cost 16.
The following bodies use the macros above. Each
SCORE includes both good-field guards (bound 353).

- In (a,b) guard b nonzero, form c by MUL and d by
  CUBE(b), lambda APPLY and XOR with KNOWN, then score.
  In (a,c) or (b,c), guard c nonzero and recover the
  missing systematic field by MUL(INV,c), then cube
  it and add its contribution to KNOWN, applying
  lambda only for a recovered b. These three bodies
  cost at most 16+18+7407+9271+1864+32+353=18961.
- In (a,d) form TARGET=lambda^-1*(d+KNOWN) via XOR
  and APPLY, omitting the constant map for (b,d).
  Guard zero TARGET and use ROOT with failure branch.
  For each of its three roots use COPY for the first
  or a zeta map for the others to recover the missing
  field, then MUL for c and score. The cap 152000
  covers 16+32+1864+18+120572+1+16+2*1864+3*(7407+353).
- In (c,d) guard c nonzero and form K=lambda*c^3
  using C1 for the cube. Use MUL(K,INV2,TARGET),
  APPLY of AS into ASROOT and APPLY of squaring into
  TEMP. XOR ASROOT into TEMP and EQ with TARGET,
  branching on failure. Set A0=MUL(d,ASROOT) and
  execute ROOT_POWER(A0,ROOT). If EQ(CHECK,A0)
  passes, branch to the cubic orientation. Otherwise
  APPLY of zeta^2 to A0 into TEMP and EQ with CHECK
  test the second orientation. On success APPLY of
  lambda^-e to ROOT yields BASEB. INV(BASEB,INVA)
  then MUL(c,INVA,BASEA) finds a. Jump to scoring.
  In the cubic orientation COPY ROOT to BASEA and
  similarly INV and MUL to find BASEB. At scoring
  COPY both base values into V0,V1 for the first
  score, then four zeta maps for the two rotations,
  scoring each. Both orientations preserve V2,V3.
  The longer path has at most
  16+18+11135+7407+1864+1864+32+34+7407+120539+
  34+1864+34+1864+131674+7407+1+32+4*1864+3*353
  =301741 operations. The cap is 310000.

When a root test or mandatory guard fails, branch
to inner extension. A good-field guard skips only
the failed candidate's score. Distinct arrays hold
each retained operand and the named macro scratch.

## Counted split chart bodies

Each split chart outer visit copies KEY0 to v's V
array, branching past a failed good-field guard.
The fixed chart indices select entries t0,t1 from
the matrix above, using at most one base XOR and
one XI_MAP. Preserve t0 in RATIO and compute
S=t1^-1 in INVCOEF by the fixed binary exponent
2^250-2 using fewer than 500 BMUL. Replace RATIO
by S*t0. For the root read the two known chunks
p0,q0 from Y, compute the missing coefficient
RATIO*p0+S*q0 by two BMUL and XOR, and assemble u.
Compute c by MUL and d from the two CUBE results,
lambda APPLY and XOR. This initializes all V arrays.

For each direction i=0..499, let f=X^(i mod 250).
For i<250 the known u coefficient changes by f
and the other by RATIO*f. For i>=250 only the
unknown coefficient changes, by S*f. Use one
BMUL and literal placement into DELTA for the
full change e and one MUL for v*e. Store e,v*e
in DIRECTIONS+16*i. The monomial f is built by
selecting its 64-bit word via shift and mask
of the index, and forming its bit by shl(1,off).
Use a runtime direction counter with bounded
test, branch, increment and jump.

For the matching affine slot set
F1=alpha*e^2, F2=alpha*e, F3=alpha*e^3 by square
and cube work and lambda APPLY where needed.
In the first half copy F1,F2 to LIN,QUAD.
In the second half use fixed theta,theta^2 maps
respectively. For each half process 32 bytes
with a bounded loop. At the first byte overall
entry zero is F3; in all other bytes it is zero.
For each literal bit r take the column LIN XOR
QUAD, using zero for the last six padded positions
in each half. Two loads, XOR and optional select
supply each output word of the column. For each
literal k<2^r and word offset, add the old row
offset, load, XOR the column, add the new offset
and store into entry k+2^r. Advance LIN by X and
QUAD by X^2 componentwise by BSHIFT_1 on the
two halves after each bit. This gives the
columns alpha*(e*beta^2+e^2*beta) for the basis
elements beta in order, plus the constant once.
The slot equals derivative_entry and evaluates
alpha*((u+e)^3+u^3). Advance the running
byte base by 2048, giving the specified fixed
layout in AFFINE for every direction.

Each byte uses up to 8*(32+90+32) for column
setup and shifts, 255*8*5 for entry expansion,
and at most 100 for the zero row, addressing,
loop comparisons, increments and jumps.
The total is below 14000. Outside the byte
loops, each direction uses one BMUL, two
field products for v*e and e^3, a square and
the constant maps, plus row stores, addressing,
and the direction loop. This fits below 100000
per direction. The coefficient inverse and
root prep fit below 1500000 more operations.
Thus the entire outer setup has bound
500*(64*14000+100000)+1500000<650000000.

At an inner visit the V arrays already match the
candidate. Guard and score for at most 353.
Before FLIP on descent save u,c,d to the row
SAVED+24*depth of the new depth by mul, add,
and add, load and store per literal word (74).
Use the slot AFFINE+i*AFFINE_STRIDE to update
d from the old u. As in APPLY, load each old
u byte and accumulate its row into the eight
old d words. Use mul/add for the slot base,
eight loads for the initial d, and an extra
add per byte for its slice offset. This costs
2+8+1864+64=1938. XOR both fields in
DIRECTIONS+16*i into u,c at cost 2+16*5=82
(mul/add for the row base, then add, two
loads, XOR and store per word). Then FLIP
advances KEY1. On ascent restore the three
fields from the current save row before
decrementing depth, with the same 74-operation
copy and the subsequent undo FLIP. The invariant
holds by root initialization and the derived
linear and cubic change identities. Counting
visit, descent and return yields 353+74+1938+
82+74=2521 <=3000 per inner node, apart from
the key traversal. Every prefix of an early
exit costs no more.

The field macros use separate scratch from
retained fields, directions and tables.
The program needs at most sixteen field
registers and twenty reused arithmetic or
address registers per active macro. Add two
depths, inner cap, BestD and table counters,
which stay within 64 registers.

## Traversal and unconditional cost

Unpack the 32 input pieces of the four fields
into Y0..Y3. Each piece loads its first word
and shifts right by its literal offset; on
crossing, load the next word, shift left and
OR. Mask when shorter than 64 and store to
its padded slot. This takes at most seven per
piece or 224 total. Move 2001 to BestD.

For every literal chart initialize KEY0
from its received outer field, outer depth
to zero and NEXT[0] to zero (17). Each
successful outer visit initializes KEY1 from
its received full field or the two selected
received halves, resets inner depth, computes
inner cap as chart cap minus outer depth and
stores zero to INEXT[0] (18). The root core
fields of a split chart already used these
same received halves in its setup.

For either tree use its KEY, depth, cap, NEXT,
LAST and labels. The inner DONE jumps by its
guard directly to the outer EXTEND, and the
outer DONE is the next chart or final output.

    VISIT: perform the stated inlined visit
    EXTEND:
        g := eq(depth,cap)
        branch if g to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        g := eq(i,500)
        branch if g to BACK
        j := add(i,1)
        store(adr,j)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        [split inner updates before the flip]
        FLIP(i)
        jump to VISIT
    BACK:
        g := eq(depth,0)
        branch if g to DONE
        adr := add(LAST,depth)
        i := load(adr)
        [split inner restoration before lowering depth]
        depth := sub(depth,1)
        FLIP(i)
        jump to EXTEND

In FLIP compute half=div(i,250), pos=mod(i,250),
baseword=shl(half,2), word=add(baseword,shr(pos,6)),
off=and(pos,63), mask=shl(1,off), and address
add(KEY,word). Load, XOR the mask, and store.
This is eleven primitives including the shift
inside word. Indices are 0..499 and padding
is untouched. A successful extension costs
25, a failed one at most six, each nonroot
backtrack 17 and root backtrack two. Summed
this is within 50 per node. Incrementing the
chosen index and storing it in NEXT before
descent visits exactly the increasing subsets
up to the cap. LAST records the flip for undo.
For radius r use O_r=sum_(i<=r) binom(500,i)
outer visits and across all outer visits
N_r=sum_(i<=r) binom(1000,i) inner visits,
because the sum of the two mask sizes is
at most r. Early outer guards reduce counts.

Pack the best fields' sixteen pieces into
their original systematic intervals low first.
Load each piece, store directly if aligned,
otherwise shift left, load, OR and store.
If it crosses also shift the original piece
right by the complementary offset and load,
OR and store to the next word. The disjoint
increasing intervals and initially zero
output guarantee unchanged earlier bits and
zero final padding. This writes all sixteen
output words with cost at most 144.

Fixed work is bounded by 224+14*17+144<1000.
Per chart the bound is (C+50)*N_r+650000068*O_r,
where C is 20000,20000,152000,20000,152000,
310000 for the six ordinary charts at r=105,
or 3000 for each of the eight split charts
at r=106. The outer coefficient includes
650000000 work, 18 to initialize inner
search and 50 control. Their sum plus 1000
gives W in the manifest. The declared run
checks W^250<=2^124993, field and memory
conditions, the exact error result and
representative chart and direction identities.
