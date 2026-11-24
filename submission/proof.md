# Fourteen information charts in a quadratic field tower

## Encoder and eligible messages

Let F be the field of order 2^250 in the polynomial basis
over the binary polynomial

    P=X^250+X^59+X^40+X^32+1.

The declared calculation checks X^(2^250)=X modulo P and
both gcd(X^(2^125)+X,P)=1 and gcd(X^(2^50)+X,P)=1.
For an irreducible factor of degree j, Frobenius has
order j on its quotient field. Its jth iterate fixes
every element by the order of the multiplicative group,
and a smaller positive iterate cannot do so by the
polynomial root bound. If an iterate fixes the residue
of X, it fixes every element since X generates that
quotient. Thus the first check forces all factor degrees
to divide 250, and the gcd checks exclude every proper
divisor. This proves P irreducible.

Form K=F[theta]/(theta^2+theta+xi) with xi=X^191.
The absolute binary trace of xi in F is one, as verified
by summing its 250 successive squares including xi.
Every u^2+u has absolute trace zero by telescoping,
so the quadratic has no root in F and K is a field of
order q=2^500. Represent each K element v=v0+theta*v1
by the 250 coefficient bits of v0 followed by the 250
coefficient bits of v1. Let lambda=X, represented by
integer 2 in F and in K. For h=(q-1)/3 put zeta=lambda^h.
The declared calculation verifies zeta!=1 and zeta^3=1.
In particular lambda is not a cube.

Encode the 1000-bit message as two K elements a,b in
sequence. The four consecutive 500-bit transmitted
groups are

    (a,b,c,d),   c=a*b,   d=a^3+lambda*b^3.

The systematic coordinates make the encoder injective
on all messages. Call an element v good when v0,v1 and
v0+v1 are all nonzero. Restrict candidate messages to
those with a and b good. Among the eligible candidates,
choose one minimizing the Hamming distance to the
received 2000-bit word, breaking ties by smallest
integer a+2^500*b in the given binary representation.
Eligibility is the union of these chart conditions.

- For each of the six pairs of full K groups, their total
  Hamming distance to the received groups is at most 105.
- Choose one full systematic group (a or b), either one
  of the two F halves of the other systematic group,
  and either one of the two F halves of c. Their total
  distance is at most 106. These are eight charts with
  1000 binary positions each.

The six full pairs also have exactly 1000 positions each.
There is always an eligible candidate. In each systematic
K group, if both halves are zero, change them to 1 and X
with two flips. If a single half is zero, flip a bit there
to 1 or X distinct from the nonzero half. If both are equal
and nonzero, flip one bit of one of them choosing a mask 1
or X that does not make its half zero. A good group already
requires no changes. Thus at most four flips make both
systematic groups good, meeting their chart cap.

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

## Admitted noise volume and error bound

Set Q=2^1000 and G=((2^250-1)*(2^250-2))^2. G is the number
of good messages because for each good field value one
may choose its first half nonzero and its second half
nonzero and different. For any t<=240, let E_t be the
set of binary errors of total weight at most t that meet
at least one of the stated chart caps. Write V_t=|E_t|.
This set depends only on the noise vector. For a received
word y let m_y be the number of good candidate codewords
within distance t, equivalently the number of good messages
whose error to y is in E_t. Then sum_y m_y=G*V_t.

Here is the exact calculation of V_t. Denote the full
group noise weights by A,B,C,D, and the smaller of the
two half weights in any group X by v_X. Let

    s=A+B,  p=C+D,  m=min(A,B),
    u=min(A+v_B,B+v_A).

All charts fail if and only if s>105, p>105,
m+C>105, m+D>105, and u+v_C>106. In particular, for total
s+p<=240 it suffices to count failures with s from 106
through 134 and p from 106 through 240-s.

For a group of total weight w, the number of masks with
both halves at least k equals

    T(w,k)=sum_(j=k..w-k) binom(250,j)*binom(250,w-j),

using T(w,k)=binom(500,w) for k<=0 and zero for k>w/2.
The relevant w values are at most 240, so the binomials
are all in range. A prefix sum over j gives these counts
exactly. For fixed A<=B the number of ordered masks with
u>=r equals T(A,r-B)*T(B,r-A). Differences between r and
r+1 give equality counts. The range of u extends from A
to min(A+floor(B/2),B+floor(A/2)); the calculation doubles
counts when A<B to account for interchanging the groups.
Fixing A,B,u and total parity weight p, the count failing
both mixed caps and all split caps is

    sum_C T(C,107-u)*binom(500,p-C),

over C>=0 with C>105-m and p-C>=0 with p-C>105-m.
Together with s>105 and p>105 this counts exactly all
failure masks. Subtract each shell's failure count from
binom(2000,s+p) and cumulate to obtain V_t.

To bound collisions, for each error in E_t allow the
competitor to have any noise of total weight at most t.
This only enlarges the number of ordered distinct pairs.
For a fixed first error and fixed systematic competitor
noise of weight j, the parity noise lies in an ordinary
1000-bit ball of radius ell=t-j. Let J be the volume of
this ball. Write D0=42. A binary parity mask is balanced
for this bound when at least one of its two 500-bit parts
has at least D0 ones and D0 zeros. Such a mask has at
least D0 ones and zeros across the parity ball. The
signed shell count of its mask of weight w is the
coefficient of z^k in (1-z)^w(1+z)^(1000-w).
Pair D0 factors of each sign. The resulting signed
coefficients are bounded in absolute value by the
coefficients of (1+z^2)^D0*(1+z)^(1000-2*D0) by
the triangle inequality. Thus the absolute signed
ball count is at most H, the sum of those coefficient
bounds through ell. Translation by the fixed first
parity error multiplies the signed count by a sign.
The affine hyperplane selected by the direction
meets the translated ball in at most (J+H)/2 values.
For an unbalanced mask use the trivial bound J.
The number of such directions is at most

    Kbad=15*(2*sum_(i=0..D0-1) binom(500,i))^2

by the multiplicity proof in either chosen binary basis.

The systematic noise of a competitor determines the
direction when the first error is fixed, by binary XOR.
For each weight class j there are N=binom(1000,j) such
patterns, with at most min(N,Kbad) unbalanced directions.
The zero direction is excluded but can increase the
upper count if included in N. From the mass 2/Q on the
affine hyperplane, the expected number of competitors
over all uniform first messages is at most A_t/Q,
where A_t sums N*(J+H)+min(N,Kbad)*(J-H) over j=0..t.
The signed coefficient bounds satisfy H<=J because
their polynomial replaces (1+z)^(2*D0) factors by
(1+z^2)^D0, which is coefficientwise no larger.
Restricting either message to the good set can only
decrease the count. Thus sum_y m_y^2<=V_t*(G+A_t).

Let F_t=|{y:m_y>0}|/Q. By Cauchy-Schwarz,
F_t>=V_t*G^2/(Q*(G+A_t)). For integer l>=1, the inequality
(x-l)(x-l-1)>=0 for positive integer x gives
l*(l+1)>=(2*l+1)*x-x^2. It also holds with zero
left side at x=0. Summing this form over y and using
both moment estimates yields
F_t>=V_t*(2*l*G-A_t)/(Q*l*(l+1)).
The exact calculation uses l=1+floor(A_t/G), floors
both bounds, and takes their maximum with the previous
lower bound. Let D_t be this nondecreasing sequence,
with D_(-1)=0.

The normalized cumulative number of received words with
winning distance at most t is precisely F_t. Each
received word with winning distance t contributes
Q^-1*9^(2000-t)/10^2000 to average success, because its
selected message has exactly one corresponding noise
pattern including on tied distances. For t=0..240 the
success lower bound is the sum of (D_t-D_(t-1))*9^(2000-t)
divided by 10^2000. Summation by parts gives nonnegative
coefficients for each cumulative D_t, so the actual
cumulative counts can only increase this expression.
The certificate finds one minus this lower bound at most
11166727434631111/200000000000000000. This normalization
includes the entire uniform message distribution.

## Solving the six ordinary charts

The chart for a specified pair of group values produces
exactly its compatible (a,b) with a,b nonzero. Chart pairs
are visited in lexicographic pair order. On each chart
enumerate all two-group received perturbations of weight
at most that pair's cap. Apply the good-message test to
each resulting candidate. Duplicate messages do not alter
the minimum or numerical tie rule.

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

## Solving the eight split charts

Take v to be the known full systematic field, the
other field u to be missing one half, and c=v*u.
Multiplication by v in the two-dimensional F basis is

    [v0    xi*v1]
    [v1 v0 + v1].

This follows from theta^2=theta+xi. Every entry of this
matrix is nonzero when v is good. Let g be the index
of the known half of u and h0 the index of the known
half of c. Denote their known values by p0,q0. In row
h0 of the matrix let t0 be column g and t1 the other
column. Then the unique missing half is
S*(q0+t0*p0), where S=t1^-1. This gives an explicit
candidate u, which is retained exactly when u is good.
The last parity field is Dbase+alpha*u^3, where
Dbase=v^3 and alpha=lambda when the known full group
is a, or Dbase=lambda*v^3 and alpha=1 otherwise.
The union of these split charts and the ordinary pair
charts, with the good-message check, is precisely the
eligible set of the decoder. Indeed every admissible
message supplies its own known field bits as one
enumerated key, and each recovered message matches
all bits of its enumerated key.

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
160 primitives. For each of the sixteen a,b words
in low-to-high order, load the corresponding BEST
word, compare by ltu and eq, and update the smaller flag
using select whenever words differ. This costs 64
and implements the numerical tie rule because the word
order respects the coefficient order. Use ltu, eq, AND,
OR and select for the comparison and update of BestD.
For each systematic word load candidate and prior best,
select and store to BEST. This costs 64 more, giving
293. Good tests and SCORE together cost at most 353.

## Ordinary chart bodies

The ordinary chart uses two 500-bit keys. The outer
group is the first named except in (c,d), whose outer
group is d. At the outer visit copy KEY0 to its field
among V0..V3, guard zero and, for the six charts in
pair order, precompute the following.

- For (a,b) and (a,d) compute a^3 into KNOWN.
- For (a,c) compute a^-1 into INV and a^3 into KNOWN.
- For (b,c) compute b^-1 into INV and lambda*b^3 into KNOWN.
- For (b,d) compute lambda*b^3 into KNOWN.
- For (c,d) compute d^-1 into INV and its square into INV2.

These operations cost less than 200000 per outer
visit, including copy and guard. At an inner visit
copy KEY1 to its field among V0..V3 at cost 16.
Then do the following. Each stated SCORE includes the
two good filters, and a filter that fails skips just its
own SCORE before the next candidate or inner extension.

- For (a,b) guard b nonzero, form c by MUL(a,b)
  and d by CUBE(b), lambda APPLY and XOR with KNOWN.
  Score (a,b). For (a,c) and (b,c) guard c nonzero,
  recover the missing systematic value by MUL(INV,c),
  and form d by CUBE, applying lambda only for a
  recovered b, then XOR with KNOWN. Score the result.
  Each of these layouts costs at most
  16+18+7407+9271+1864+32+353=18961, below 20000.

- For (a,d), form TARGET by XOR(d,KNOWN) followed
  by the lambda^-1 map in place. For (b,d), omit that
  map. Guard zero TARGET. Use ROOT(TARGET,ROOT),
  with an additional failure branch. For the base
  root COPY into the missing systematic field, and for
  each of the two rotations use APPLY of zeta or
  zeta^2 to ROOT into the missing field. Form c by
  MUL(a,b) and score for each root. This uses at most
  16+32+1864+18+120572+1+16+2*1864+3*(7407+353)
  =149527, below 152000.

- For (c,d), guard c nonzero and form K=lambda*CUBE(c).
  Use C1 as the cube result before applying lambda
  into K. Use MUL(K,INV2,TARGET), APPLY of AS to
  TARGET into ASROOT, and APPLY of squaring to ASROOT
  into TEMP. XOR ASROOT into TEMP and use EQ(TEMP,TARGET)
  and a failure branch. MUL(d,ASROOT,A0) forms the
  candidate A. ROOT_POWER(A0,ROOT) puts its root test
  in CHECK. Use EQ(CHECK,A0) and branch to ordinary
  orientation if true. Otherwise APPLY zeta^2 to A0
  into TEMP, then EQ(CHECK,TEMP) and branch on failure.
  On success APPLY of lambda^-e to ROOT yields BASEB.
  Use INV(BASEB,INVA) and MUL(c,INVA,BASEA) and jump
  to shared scoring code. The ordinary orientation
  instead uses COPY(ROOT,BASEA), INV(BASEA,INVA) and
  MUL(c,INVA,BASEB), then falls through to scoring.
  In the shared code copy BASEA and BASEB to V0 and V1
  for the first score. For each rotation use two zeta
  maps of BASEA and BASEB before scoring. The c and d
  arrays are unchanged. The more expensive orientation
  uses at most

      16+18+11135+7407+1864+1864+32+34+7407+120539
      +34+1864+34+1864+131674+7407+1+32+4*1864+3*353
      =301741

  primitives, below 310000. Failed guards only reduce
  the count, and strict improvements and ties use the
  same bound. All guards, jumps, and field operations in
  these sequences were charged.

## Split chart setup and table generation

For a split chart, KEY0 is the full known field v.
KEY1 consists of p0 and q0 in that order as two
four-word halves. The relevant received halves are
copied into KEY1 at inner initialization. At the
outer visit copy KEY0 into its appropriate V0 or V1
field and guard good, jumping to the outer extension
if it fails. Find t0,t1 by the fixed matrix row and
column choice. At most one XI_MAP and one base XOR
form the entries. Copy t0 to RATIO and t1 to INVCOEF
using first halves of those arrays. Compute S=t1^-1
there by a fixed binary square-and-multiply chain for
the exponent 2^250-2, using BMUL for squaring and
products and a preserved copy of t1, disjoint from
MULT. This takes fewer than 500 BMUL plus copies.
Also compute t0*S in the first half of RATIO. The
split initialization can now use RATIO*p0+S*q0
to recover the missing half of u at the root. Read
p0,q0 from the original received groups. Copy p0
to the known half of u, form c by MUL and d using
one CUBE for v, one CUBE for u, one lambda map on
the appropriate cube and XOR. Store the four groups
in V0..V3 at this point, including u and c.

Compute the direction tables for this outer key
before the inner traversal. For i from 0 to 499 let
f=X^(i mod 250) in F. The direction uses a flip of
the known systematic half p0 when i<250, and of q0
otherwise. The linear change e in the full field
u has its gth half f in the first case and zero in
the second. Its other half is respectively (S*t0)*f
or S*f. Use one BMUL, the four-word representation
of f, and literal copies into DELTA to form e.
The direction row stores e and v*e, generated by a
MUL. To obtain f, use mod, shift and AND to find
the word and bit of i mod 250, then one shl from
literal 1 and per literal output word an equality,
select and store. The simple branch i<250 chooses
which literal version to execute. Each row is at
DIRECTIONS+16*i and uses sixteen dynamic stores
with address increments or adds.

For the matching affine table, form F1=alpha*e^2,
F2=alpha*e, F3=alpha*e^3. Use APPLY for squaring,
one additional MUL for the cube or CUBE and use
the lambda map on each expression for alpha=lambda.
In the first half start LIN=F1, QUAD=F2. In the
second half start LIN=theta*F1, QUAD=theta^2*F2
by the fixed maps. Process 32 bytes for each half.
The loop over half choice can be emitted twice with
literal initial maps and offsets. For each byte,
first store entry zero as F3 in the first slice
of the whole array and zero otherwise. For each of
its eight literal bits take the column LIN XOR QUAD,
except that for the last six bits of each half use
the zero column. Obtain each output word by two
loads and XOR, optionally select zero using the
byte comparison. For each literal k from 0 through
2^r-1 for bit index r and each literal word, add its
offset to the slice base, load entry k, XOR the
column word, add the destination offset, then store
to entry k+2^r. Keep each reused column word in a
register while emitting the literal k copies.
This expansion establishes all 256 entries by
induction over r. Advance LIN by X and QUAD by X^2
componentwise, with BSHIFT_1 once or twice on each
of their two halves after each bit, loading and
storing the corresponding words. The direction
table is exactly the map

    u -> alpha*(e*u^2+e^2*u+e^3),

since each basis bit contributes alpha*(e*beta^2+e^2*beta)
at its basis coordinate beta and the constant is in
the first slice only. This is the table determined by
derivative_entry. Use counters with fixed bounds 500
for directions and 32 for each half's byte loop.
Emit literal comparisons with the bounds, branches,
increments and jumps. Advance a slice cursor by
2048 per byte from its start at AFFINE, giving
the index AFFINE+i*AFFINE_STRIDE for direction i.

The column step uses at most eight groups of
8*4+6*15+32 operations each, including LIN and
QUAD loads and stores. The entry expansion uses
255*8*5. An allowance of 100 more covers row zero,
padding and loop comparisons, address maintenance
and branches. This is less than 14000 per slice.
The remaining work per direction is smaller than
100000. For example it uses one BMUL, one MUL for
v*e, one square and one CUBE for the seeds, and up
to five constant APPLY maps. All row addressing,
copies, selection and initialization also fit below
this allowance. Thus the direction loops cost at
most 500*(64*14000+100000) plus fixed initialization.
The coefficient inverse and root initialization
use fewer than 1500000 further operations. The
uniform cap 650000000 covers the entire outer visit,
including its guard and all table writes.

At an inner split visit all four V arrays already
hold the candidate. The visit performs the two good
tests and score, of cost at most 353. Before FLIP on
descent, save u,c,d into SAVED+24*depth of the new
depth. Emit mul, add, then one add, load and store
per literal word, costing 74. Use mul and add to get
the affine table base for direction i and apply it
to the old u, with accumulators initialized by
loading the eight old d words. Use APPLY with one
extra add per byte to add the literal slice offset
beyond the direction base. Store the updated d.
This costs 2+8+1864+64=1938. It is performed before
changing u. Read the row of DIRECTIONS and XOR its
two fields into u,c. The cost is 2+16*5=82 including
the mul and add for the row and add, two loads,
XOR and store per literal word. Thus u and c reflect
the changed p0,q0, and d is its correct cubic parity.
Then execute the usual KEY1 flip. On ascent restore
u,c,d from the current depth's save row before
decrementing depth, with the same 74-operation
addressing and copy. Flip the indicated key bit back.
The invariant holds at the root by initialization
and at every descendant or return by these updates.
Together, the visit and one entry and exit cost at
most 353+74+1938+82+74=2521, below 3000.
Traversal control and key flips are charged separately.

The named arithmetic scratch arrays used inside MUL
are distinct from the retained chart state. The AS,
CUBE and geometric macros use the temporary arrays
specified in their descriptions. In split setup DELTA,
RATIO, INVCOEF, F1,F2,F3, LIN and QUAD do not alias
these temporaries. One field macro at a time runs,
using at most sixteen field registers and fewer than
twenty reusable arithmetic, loop and address
temporaries. The traversal retains the two depths,
the inner cap and BestD in four additional registers.
The table builder can reuse the existing scope's scratch
registers and has only its direction and byte counters
and slice cursor live across column steps. These bounds
fit within 64 registers.

## Traversal, final output and total count

At initialization extract the four received fields
into Y0..Y3 as padded word arrays. Each of the 32
pieces is a contiguous at-most-64-bit interval.
Load its first source word and shift right by the
fixed bit offset if necessary. If its live bits
cross the word boundary, load the next word, shift
left by 64 minus the offset and OR. Mask to its
width if shorter than 64, and store the piece.
Each piece uses at most seven primitives including
its store, hence total work at most 224.

Initialize the distance register BestD to 2001 by a
free move. For each chart copy its outer received
field into KEY0, move the outer depth to zero and
store zero in NEXT[0], costing 17. At each outer
visit run the chart's outer body and, if its guard
passes, initialize the inner traversal. For ordinary
keys copy the inner received field into KEY1. For
split keys copy the two received halves into the
two halves of KEY1. Either operation costs 16.
Initialize the inner depth to zero by a move,
compute the inner cap by subtracting outer depth
from chart cap, and store zero at INEXT[0]. This
costs 18 in total. For split charts the outer
body initializes its V arrays using those same
received halves, so KEY1 initialization does not
change the candidate fields.

For each outer and inner traversal the following
template is emitted. The pair of depth registers is
distinct. In an outer copy use arrays KEY0,NEXT,LAST
and the fixed chart cap. In an inner copy use
KEY1,INEXT,ILAST and the computed inner cap. The
inner DONE is the outer EXTEND, and outer DONE
falls through to the next chart or final output.

    VISIT:
        execute the corresponding inlined action
    EXTEND:
        g := eq(depth,cap)
        branch if g to BACK
        addr := add(NEXT,depth)
        i := load(addr)
        g := eq(i,500)
        branch if g to BACK
        j := add(i,1)
        store(addr,j)
        depth := add(depth,1)
        addr := add(NEXT,depth)
        store(addr,j)
        addr := add(LAST,depth)
        store(addr,i)
        [in split inner copy, save and update here]
        FLIP(i)
        jump to VISIT
    BACK:
        g := eq(depth,0)
        branch if g to DONE
        addr := add(LAST,depth)
        i := load(addr)
        [in split inner copy, restore here]
        depth := sub(depth,1)
        FLIP(i)
        jump to EXTEND

A branch not taken falls through. FLIP emits

    half := div(i,250)
    pos := mod(i,250)
    baseword := shl(half,2)
    word := shr(pos,6)
    word := add(baseword,word)
    off := and(pos,63)
    mask := shl(1,off)
    addr := add(KEY,word)
    v := load(addr)
    v := xor(v,mask)
    store(addr,v)

where div and mod are unsigned integer operations.
The 11 primitives never affect padding. Each
successful extension costs 25 aside from the
charged split updates, each nonroot BACK costs
17, a failed extension costs at most six, and
root BACK costs two. Hence the traversal overhead
is at most 50 per node summed over either tree.
Each visited mask uses increasing indices, since the
parent NEXT is advanced and the child NEXT starts
at the next greater bit. The recording in LAST and
the symmetric flip on return recover the parent
mask. By induction the traversal visits exactly the
subsets of size at most its cap. The outer and
inner masks run over four F halves in either chart
type, with exactly 500 bits assigned to each key.
For radius r the outer node count is
O_r=sum_(i<=r) binom(500,i). Across these outer nodes,
the inner cap accounts for outer weight and the
total potential inner count is N_r=sum_(i<=r) binom(1000,i).
A failed outer guard can only decrease that count.

At final output, read the eight pieces of the best
a and likewise the eight pieces of the best b. Each
piece is packed into its original 250-bit-half layout,
with the next half immediately following the previous
one. If the destination is word aligned, emit a direct
store after loading. Otherwise shift the piece left by
the bit offset, load the output word, OR and store.
When live bits cross a boundary, also shift the
original piece right by 64 minus the offset and load,
OR and store to the next word. The pieces cover the
1000 output bits with disjoint intervals, and their
increasing order means aligned stores can target only
words with no earlier bits. The unused final 24 bits
remain zero. This uses at most nine operations per
piece, for at most 144. Every output word is written.

Fixed work is at most 224+14*17+144=606, below 1000.
For a chart of cap r and per-inner bound C combine
(C+50)*N_r and 650000068*O_r, which includes the
outer precomputation, inner initialization and
outer traversal overhead. The six ordinary bounds
are (20000,20000,152000,20000,152000,310000) at
radius 105; each split bound is 3000 at radius 106.
Using the integer binomial sums and adding 1000
gives the stated integer W in submission.json. The
declared computation checks W^250<=2^124993 exactly.
It also checks the memory extent, field identities,
admitted volume, ordinary and split chart equations
and sample direction-update identities. The proofs
of the fully expanded program and the exact count
apply without executing the full search.
