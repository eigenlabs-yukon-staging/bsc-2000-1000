# Eight split information charts in a quadratic tower

## Code and decoder

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

Send a,b,c,d=(a,b,ab,a^3+lambda*b^3) in four
successive groups, with the first two groups
carrying all 1000 input bits. A field value is
good when its two F coefficients are nonzero
and distinct. Search the eight choices of
one full systematic group, either half of the
other systematic group, and either half of c.
Keep good a,b from each known-coordinate key
whose combined distance is at most 103.
Additionally consider the fixed fallback
a=b=1+theta*X. It is good and guarantees a
nonempty set. Choose minimum full distance,
breaking ties by numeric a+q*b.

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

## Admitted noise and success

Let V_t count errors of total weight at most t
that meet at least one split cap, excluding the
fallback from this volume calculation. For total
weights A,B,C,D in the four groups, and v_X the
smaller half weight of group X, put
u=min(A+v_B,B+v_A). The admission condition
is u+v_C<=103, exactly the minimum over all
eight charts.

For a group of total w, the count with half
weights j,w-j is binom(250,j)*binom(250,w-j).
Group these counts by their minimum and take
suffixes T(w,k) for minimum at least k. At
fixed A<=B the count with u>=l is
T(A,l-B)*T(B,l-A). Differences for consecutive
l yield equality counts for u, and double
when A<B. The range A..min(A+floor(B/2),
B+floor(A/2)) contains all u. Sum over A+B=s
to get its exact distribution per systematic
total s. For parity total p, combine the
minimum-half counts of C with binom(500,D) over
all C+D=p. Take prefix sums in v_C. For any
total weight t, combine the exact u bins of
systematic total s and the parity prefixes
through 103-u for p=t-s. Summing over s and
cumulating over t gives V_t for every t<=240.

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
success. The winning-distance cumulative count divided
by Q is at least F_t, since the fallback can
only lower the selected distance. Hence
sum_(t<=240)(D_t-D_(t-1))*9^(2000-t)/10^2000
lower-bounds success, by summation by parts
with nonnegative coefficients. The declared
computation obtains error at most
109323247683769923/1000000000000000000.

## Split chart recovery

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
words and a K element takes eight. The 30 arrays in ARRAYS start at 48 with
eight words each. The four traversal stacks
NEXT,LAST,INEXT,ILAST have 104 words each at
288,392,496,600. The 16-row four-word MULT
scratch starts at 704. The fixed two- and
sixteen-entry reduction arrays start at 768
and 770. XI_MAP at 786 has 32 slices of 256
four-word entries for xi multiplication.
MAP_BASE at 33554 has four maps of 64 byte
slices and 256 eight-word entries per slice,
covering squaring and multiplication by
lambda,theta,theta^2. There are 32 bytes in
each half, with two live positions in the
last byte. xi_entry, map_columns and map_words
give every entry, assigning zero to padding
columns. DIRECTIONS at 557842 has 500
sixteen-word rows, SAVED at 565842 has 104
24-word rows, and AFFINE at 568338 has 500
direction slots each the size of one fixed
K map. All scratch and read-only ranges are
disjoint and end at exclusive address
66104338, below 2^32.

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

CUBE uses the square map into X2 and MUL with
the original operand, at 9271 operations.
Ensure those operands are distinct from X2.
The root and table seed work uses the named
other arrays. Base inversion below uses a
preserved F word-pair and BMUL.

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

## Counted split chart bodies

Each split chart outer visit copies KEY0 to v's V
array, branching past a failed good-field guard.
The fixed chart indices select entries t0,t1 from
the matrix above, using at most one base XOR and
one XI_MAP. Preserve t0 in RATIO and save t1 in BASEPIVOT and compute
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
82+74=2521 <=2600 per inner node, apart from
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
piece or 224 total. Move 2001 to BestD. Set the four V arrays
to all thirty-two literal padded words of
codeword(FALLBACK,FALLBACK) by stores and
score the candidate before any charts.

For every literal chart initialize KEY0
from its received outer field, outer depth
to zero and NEXT[0] to zero (17). Each
successful outer visit initializes KEY1 from
the two selected received halves, resets
inner depth, computes
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

The fixed input, fallback and output work is at
most 224+32+353+144<1000. Count 20 for each chart's
start and (2600+50)*N_r+650000068*O_r at r=103.
The outer term covers its table and arithmetic work,
18 for inner initialization and 50 traversal.
Summing for the eight charts and adding 1000
gives W in the manifest. The certificate checks
W^100<=2^48839, the exact error bound, field and memory
conditions and sample direction updates. The full
decision and operation proof applies on every input
without executing the full enumeration.
