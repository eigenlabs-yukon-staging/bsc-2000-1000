# Incremental two-sided field-cube decoding

## Field code and candidate rule

Let m=1000 and Q=2^m. Use K=F2[X]/(P) with
P=X^1000+X^53+X^45+X^34+1. In the polynomial basis, coefficient
of X^i is bit i. The declared computation verifies irreducibility
by repeated modular squaring from X. It checks X^(2^1000)=X modulo
P and gcd one of P with X^(2^500)+X and X^(2^200)+X.
An irreducible factor of degree d divides X^(2^j)+X if and only
if d divides j. In its field of size 2^d, Frobenius has order d:
its d-th power fixes the field by the order of its multiplicative
group, while any j with 0<j<d would otherwise give 2^d roots of
a polynomial X^(2^j)+X of lower degree. Fixing the residue of X
is equivalent to fixing its entire quotient field. Thus every
factor degree of P divides 1000. A proper divisor of 1000 divides
1000/2 or 1000/5, so the two gcd checks exclude all proper degrees.
It follows that P is irreducible.

Encode u as f(u)=(u,u^3). The first m coordinates are the m-bit
message and the last m are its field cube. This is an injective
encoder of 2^1000 messages in 2000 uses. For a received (a,b),
the decoder considers each v satisfying
wt(v+a)<=112 or wt(v^3+b)<=109. It minimizes the sum of the two
weights, breaking ties by smallest numerical v in the fixed
polynomial basis. The set is nonempty because it contains v=a.
The decoder is deterministic.

## A field differential identity

Fix a nonzero systematic difference d in K. For uniform u,
(u+d)^3+u^3 is distributed uniformly on a trace hyperplane.
In fact expansion gives d^3((u/d)^2+(u/d)+1).
The map z -> z^2+z is F2-linear with kernel {0,1}, and its
image has size 2^(m-1). Write Tr(z)=sum_(i=0..m-1) z^(2^i).
This is F2-linear with values in F2, since its square is itself.
It is nonzero since its defining polynomial has degree 2^(m-1)
and cannot vanish on the whole field. Also Tr(z^2+z)=0, and
Tr(1)=0 because m is even. Therefore the set of possible parity
differences is exactly
H_d={c in K: Tr(d^-3 c)=0},
each with probability 2^(1-m).

For each d the equation has a nonzero binary coefficient mask h,
whose coordinate j is Tr(d^-3 X^j). The map sending any field
element q to this trace-coordinate vector (with q in place of
d^-3) is an F2-linear bijection. Its kernel is zero: if all these
traces vanish and q is nonzero, multiplication by q is surjective,
which would make the trace identically zero. Each nonzero q has
at most three preimages under d -> d^-3, by the root bound for
the cubic equation. It follows that the number of differences
d for which wt(h) is at most D or above m-D is at most
Kbad=3*(V_D-1+V_(D-1)), where D=48 and
V_j=sum_(i=0..j) binom(m,i). The subtracted one excludes the
zero mask from the low-weight range. The high-weight range is
counted by complementing mask bits.

## Parity-ball intersections

For a parity-space binary mask h with s=wt(h)>=1 and integer
ell>=0, the signed sum of (-1)^(h dot z) over wt(z)<=ell
is the coefficient of x^ell in
(1-x)^(s-1)*(1+x)^(m-s).
To see this, the contribution at exactly degree j is the
coefficient of x^j in (1-x)^s*(1+x)^(m-s), by choosing the
positions relative to h. Cumulative coefficients divide the
generating series by 1-x.

For D+1<=s<=m-D there are at least D pairs consisting of a
factor 1-x and a factor 1+x. Pair all possible factors, say q
pairs with q>=D. Taking the absolute value coefficientwise
is bounded by coefficients of
(1+x^2)^q*(1+x)^(m-1-2q).
Since 1+x^2 is bounded coefficientwise by (1+x)^2,
this is bounded coefficientwise by
(1+x^2)^D*(1+x)^(m-1-2D).
Let C_ell be its coefficient. For the other masks the absolute
signed sum is bounded by V_ell. Define both values to be zero
at ell=-1. Centering the ball at any fixed parity vector b
multiplies its signed sum by (-1)^(h dot b), which does not
affect the absolute bound. Hence its intersection with the
trace hyperplane has cardinality at most (V_ell+C_ell)/2
for the balanced masks and at most V_ell for the other masks.
The computed C_ell is no larger than V_ell, also from the
coefficientwise bound by (1+x)^(m-1).

## Average success from occupied received words

For total weight t let E_t consist of the vectors of weight at most t
that satisfy at least one of the two conditions: systematic weight at
most R_I=112, or parity weight at most R_P=109. For a fixed received
word y, let m_y count the messages v with y+f(v) in E_t.
Let J_t be the cardinality of E_t, H_t be the number of y with
m_y>0, and F_t=H_t/Q. Every sent word has exactly J_t errors
in E_t, so sum_y m_y=Q J_t.

When m_y>0, the least admitted distance of y is at most t.
Conversely only those y have such a least distance at most t.
At a fixed y with winning distance w, the joint success
contribution, summing over uniform messages, is
Q^-1 9^(2000-w)/10^2000. It is the probability of the selected
message and its uniquely determined noise vector. Nearest
candidate ties all yield the same value. Consequently the
exact average success is
sum_t (F_t-F_(t-1)) 9^(2000-t)/10^2000, starting with F_(-1)=0.

For a fixed true error z in E_t, count competing messages
whose errors to the same received word also belong to E_t.
Let s be the systematic weight of a competitor's error.
For a fixed message difference d, this weight is wt(z_I+d),
independent of the true message u. For s<=t the allowed
parity-space ball has radius ell=t-s when s<=R_I and otherwise
radius min(t-s,R_P). It is centered at z_P for the parity
difference f(u+d)_P+f(u)_P. Write v=V_ell and c=C_ell
for the ordinary m-bit ball size and signed-sum bound from
the previous section. There are at most n_s=binom(m,s)
nonzero differences with that systematic weight. Among
them at most min(n_s,Kbad) are unbalanced. For a balanced
difference the expected indicator of membership in this ball,
averaged over u, is at most (v+c)/Q; for the other
differences it is at most 2v/Q. Thus for each fixed z the
expected count of all competing messages in E_t is at most
A_t/Q, where

    A_t = sum_(s=0..t) [n_s*(v+c)+min(n_s,Kbad)*(v-c)].

The experiment computes this integer and also
J_t=sum_(s=0..t) n_s*v, using the same ell rule.
The difference count excludes d=0 in the definition of
a competitor; using n_s as an upper bound is safe.

For each ordered pair of distinct messages counted at
a y with both errors in E_t, take the first message
as the true one and its corresponding z as the fixed error.
Then average the competitor count over that uniform true
message and sum over all J_t choices of z. This counts
each ordered pair at each y exactly once and yields
sum_y m_y(m_y-1) <= J_t A_t. It follows that
sum_y m_y^2 <= J_t(Q+A_t). Cauchy-Schwarz gives
F_t >= J_t Q/(Q+A_t).

For any integer k>=1 and positive occupancy x, use
(x-k)(x-k-1)>=0 to obtain
k(k+1)>= (2k+1)x-x^2. At x=0 the right side is zero,
so summing this inequality over occupied y and using
the collision bound gives
F_t >= J_t*(2k Q-A_t)/(Q k(k+1)).
The run takes k=1+floor(A_t/Q) and floors both rational
bounds, keeping the maximum with the previous lower
bound. Write the resulting integer lower values D_t.
They are nondecreasing and satisfy D_t<=F_t.
The success numerator is at least
sum_(t=0..240) (D_t-D_(t-1))*9^(2000-t).
Indeed summation by parts expresses the truncated sum
using positive differences of successive weights, plus
a nonnegative final weight. Substituting F_t for D_t
only increases it, and the exact success includes
nonnegative contributions after t=240. The experiment
checks exactly that one minus this lower numerator
over 10^2000 is at most
55176594175514213/1000000000000000000.

## Inverting the parity candidates

The decoder makes two passes. The systematic pass
enumerates the radius-112 ball about a and cubes each
element, considering the resulting message. The parity
pass enumerates the radius-109 ball about b and considers
all cube roots of every key z in that ball. Any
message meeting the candidate rule occurs in at least
one pass. Duplicates are allowed and do not change the
minimum or numerical tie rule.

Let h=(Q-1)/3 and e=(h+1)/3, which are integers as
checked in the declared computation. For nonzero cubic z,
z^h=1 because z is a cube of a field element and the
multiplicative group has order Q-1. Thus (z^e)^3=z.
The same equality holds at zero. For other z the computed
root fails the cube check, since equality would make it
a cube. Let Zeta=7^h, with 7 represented in the fixed
polynomial basis. The computation checks Zeta!=1 and
Zeta^3=1. For a valid nonzero root x, the three values
x, Zeta*x, Zeta^2*x are distinct roots and exhaust the
cubic by degree. At zero all three are the same correct
root. For each valid key the program scores these three
candidates against the original input.

The fixed calculation of z^e uses
e=2+112*sum_(i=0..165) 2^(6i). Compute Z2=z^2 and
Z4=z^(2^2), multiply those two, multiply by z, and
raise the result to 2^4 to form T=z^112. For
G_a=T^(sum_(i=0..a-1)2^(6i)), the relation
G_(2a)=G_a*(G_a)^(2^(6a)) doubles the length.
For a further term multiply by T^(2^(12a)).
Begin with G_1=T, and append bits (0,1,0,0,1,1,0)
to the leading one of 166. These steps compute G_166.
Multiplying it by Z2 gives the required root candidate.
The declared run checks the exact exponent identity
and the field arithmetic on samples.

## Fixed memory and primitive field macros

Represent a field value by sixteen consecutive words
with 40 live bits in its high word. Input/output memory
is as specified in the rules. Starting at 48, put sixteen-word scratch arrays in the order
YS, YP, KEY, SQ, CUBE, BEST, Z2, Z4, T, G, H, J, ROOT,
ROTATED, P1, P2. Then put two 113-word arrays NEXT and LAST
at 304 and 417. At 530 put multiplication scratch MULT, with
64 rows of sixteen words. At 1554 put the 32-word RAW scratch
array. At 1586 put HISTORY with 113 rows of 32 words.
At 5202 and 5204 put the fixed reduction tables of two and
64 words. For s=1 or 6 their entry c is the carryless product
of c with P_low=X^53+X^45+X^34+1, which fits a word.
At 5268 put the fixed directional data. At offsets
32*i through 32*i+15 are the sixteen words of X^(2i) mod P,
then at offsets 32*i+16 through 32*i+31 those of X^(3i) mod P.
There are 1000 such groups, as in direction_words().
All ranges including the input/output prefix are disjoint.

Starting at address 37268, allocate 15 fixed linear-map tables
in the order map_tags(): Frobenius exponents
1,2,4,6,12,24,30,60,120,240,246,492,498, and then multiplication
by Zeta and Zeta^2. For a fixed map F, at its table base plus
4096*j+16*d+h store word h of F(d*X^(8*j)), with d interpreted
as a byte polynomial. The map_words() generator determines all
of them by applying the indicated map to basis coefficients
and combining linearly. The maps are linear over F2, so the
XOR of the 125 selected field table values is the full map
value. The tables end at exclusive address 7717268.
The remaining working memory is in the scratch ranges.

All loops with literal ranges in the next macros are expanded
at code generation time. Moves and constants are free. A macro
invocation is inlined with input and output base addresses
substituted as constants. The actual program's only control
primitives are written explicitly below. Words use the
arithmetic and shifts defined by bsc-ops-v1.

SHIFT_s multiplies a 16-register field value R by X^s modulo P
for literal s=1 or 6. Emit:

    c := shr(R15,40-s)
    for h=15 down to 1:
        t := shr(R[h-1],64-s)
        R[h] := shl(R[h],s)
        R[h] := or(R[h],t)
    R0 := shl(R0,s)
    R15 := and(R15,2^40-1)
    adr := add(REDUCE_s,c)
    t := load(adr)
    R0 := xor(R0,t)

The descending order preserves the original word for a carry.
The high coefficients above degree 999 are reduced using
X^1000=X^53+X^45+X^34+1. Their correction lies wholly in the
low word because s<=6. This costs 1+45+1+1+3=51.

APPLY(F,A,D) moves sixteen accumulator registers R to zero.
For each byte index j=0..124, load from A+floor(j/8), then
shift right by 8*(j mod 8), AND with 255, shift left by four,
and add the literal table base for F plus 4096*j. For each
h=0..15 add h to that address, load its word, and XOR it into
R[h]. Store R to D with sixteen literal stores after the last
byte. This costs 125*(5+48)+16=6641. Output writes occur only
after all input reads, so A=D is valid.

MUL(A,B,D) loads A into field registers. Store zero and this
value to rows 0 and 1 of MULT, costing 16+32 for the loads and
stores. Apply SHIFT_1 five times in succession and store the
results in rows 2,4,8,16,32. For each other nonzero row d in
increasing order write the XOR of the row for its least set bit
and the row for the remaining bits. This requires two literal
loads, one XOR and one literal store for each of its words.
Thus row d is A times polynomial d of degree at most five.

Move R to zero. For each literal six-bit digit index j from 166
down to zero, let (w,s)=divmod(6j,64) at code generation.
Load word B+w and shift right by s. Where s>58 also load
B+w+1, shift it left by 64-s, and OR it into the result,
at three additional operations. AND with 63. Shift the digit
left by four and add MULT to form the row address. Use SHIFT_6
to multiply R by X^6 and for each word add its literal index
to the row address, load and XOR into R. This is Horner
evaluation of the field product. The highest digit uses zero
padding above bit 999. At the end store sixteen result words
in D. The 167 digits include ten crossing cases, exactly
those where 6j mod 64 is 60 or 62. The resulting primitive
count is 16+32+5*(51+16)+57*64+167*(5+51+48)+3*10+16=21445.
The macro reads A and B before the final result stores,
so D may alias either. MULT is separate from all operands.

## Incremental cubing during systematic search

The systematic pass keeps SQ and CUBE equal to the square and
cube of KEY. On initialization, after copying YS to KEY,
emit APPLY(1,KEY,SQ) and MUL(SQ,KEY,CUBE). This costs once
6641+21445=28086. If the next flip is coordinate i, let
e=X^i. In characteristic two, the new values are

    new_sq = SQ + e^2
    new_cube = CUBE + e*SQ + e^2*KEY + e^3.

All right-hand sides refer to the current parent. We compute the
three products by powers of X with the macro MONOMIAL(A,i,D)
below, along with the fixed e^2 and e^3 directional data.
This uses polynomial shifts and reduction by the sparse P.

MONOMIAL computes X^i*A for a literal input base A and output
base D, and a register i in 0..999. Clear all 32 RAW words by
literal stores. Set q=shr(i,6), s=and(i,63) and v=sub(64,s).
For each literal j=0..15 emit the following instructions.

    u := load(A+j)
    lo := shl(u,s)
    hi := shr(u,v)
    a := add(RAW+j,q)
    b := add(a,1)
    t := load(a)
    t := xor(t,lo)
    store(a,t)
    t := load(b)
    t := xor(t,hi)
    store(b,t)

These 32+3+16*11 primitives place all the polynomial bits of
A shifted left by i in RAW, with contributions at word
boundaries combined by XOR. A shift right by 64 yields zero,
so the s=0 case is included. RAW addresses are at most RAW+31.

Split the raw value W into its m=1000 low bits and the
higher polynomial H=W>>1000. For each literal j=0..15,
load RAW+15+j and RAW+16+j, shift right by 40 and left
by 24 respectively, and OR to get register H[j]. This uses
five primitives each. View H[-1] and H[16] as zero constants.
For each j=0..16 let T[j] start as H[j]. For each of the
three literal exponents a=34,45,53, XOR into it the two
values shl(H[j],a) and shr(H[j-1],64-a). That is twelve
primitives per j. For j<16 load and XOR RAW+j into T[j],
with an additional AND by 2^40-1 on that input at j=15.
The T registers are now the polynomial
(W mod X^m)+P_low*H, costing a further 17*12+33.

Its degree is at most 1051, since the original degree
is at most 1998 and H has degree at most 998. To fold
once more, take b=(T[15]>>40) OR (T[16]<<24), costing
three, and AND T[15] with 2^40-1. This b has at most 52 bits.
Set two registers Z0=b and Z1=0. For each of the same three
exponents a, XOR shl(b,a) into Z0 and shr(b,64-a) into Z1,
at twelve operations altogether. This forms the low and
high words of b*P_low, which has degree at most 104.
XOR Z0 and Z1 into T[0],T[1], and store T[0..15] at D.
All later coefficients vanish, so this is the required
remainder modulo P. The last stage costs 3+1+12+2+16=34.
In full MONOMIAL costs 32+3+16*11+16*5+17*12+33+34=562.
It reads A before writing D, allowing them to alias.

On a successful systematic extension, while still at the
parent depth, save SQ and CUBE in the current HISTORY row
and then update them. To save, shift depth left by five,
add the HISTORY base, then for each of the 32 words add its
literal offset, load the corresponding SQ or CUBE word,
and store at that address. The save and its inverse restore
each cost 2+32*3=98. Next apply
MONOMIAL(SQ,i,P1), MONOMIAL(KEY,i,P2), and MONOMIAL(P2,i,P2).
These provide e*SQ and e^2*KEY. Compute the directional
base by shifting i left five and adding DIRECTION_TABLE.
For each h=0..15 load SQ+h, add h to the directional
base, load that directional word, XOR and store back.
Also load CUBE+h, P1+h, and P2+h, add the literal 16+h
to the directional base, load the directional word,
XOR the four values with three XORs, and store at CUBE+h.
This costs 2+16*(5+9)=226. Only after this update flip KEY.
On a return from a systematic child and after subtracting
one from depth, restore SQ and CUBE from the parent
HISTORY row and flip KEY back. This proves the invariant
by induction on all visits and returns. The parity pass
does not use this update or HISTORY.

## Candidate evaluation and word RAM traversal

For each of the first 15 input words load and store it into
YS. For input word 15 mask its low 40 bits before storing.
These 33 primitives set YS. To initialize YP, for each
literal j=0..15 load words 15+j and 16+j, shift the first
right by 40 and the second left by 24, OR and store to YP+j.
This costs 96. Word 31 has 16 live bits, so the high word
of YP already has the required 24 padding bits zero.
Move 2001 into BestD. Before each pass, copy YS or YP
to KEY using a load and a store for each word, write zero to
NEXT[0] and move zero into depth. This costs 33 per pass.

SCORE(A,B) tests the message at array A and its parity at
array B. Start Dist and Smaller at zero by free moves.
For each literal h=0..15, emit:

    a := load(A+h)
    b := load(YS+h)
    t := xor(a,b)
    t := popcount(t)
    Dist := add(Dist,t)
    a2 := load(B+h)
    b2 := load(YP+h)
    t := xor(a2,b2)
    t := popcount(t)
    Dist := add(Dist,t)
    b := load(BEST+h)
    lt := ltu(a,b)
    same := eq(a,b)
    Smaller := select(same,Smaller,lt)

The loop's low-to-high order makes the most significant
differing word determine the numerical comparison.
Use ltu(Dist,BestD), eq(Dist,BestD), AND of the equality
with Smaller, and OR with the strict-improvement flag.
Select BestD according to the result. For each h load
A+h and BEST+h, select the former when this flag is set
else the latter, and store at BEST+h. The cost is
16*14+5+16*4=293. The first candidate has distance at
most 2000, and installs itself. Later tests implement
the exact minimum and tie rule.

In the systematic pass the VISIT action is SCORE(KEY,CUBE),
costing 293, because the incremental update already maintains
the correct cube. In the parity pass the VISIT action is
the following inlined sequence, using integer map tags for
Frobenius exponents and zeta tags as specified above:

    APPLY(1,KEY,Z2)
    APPLY(2,KEY,Z4)
    MUL(Z2,Z4,H)
    MUL(H,KEY,T)
    APPLY(4,T,T)
    copy T into G                             # 32 loads and stores in total
    a := 1                                   # literal at code generation
    for bit in the fixed chain:
        APPLY(6*a,G,H)
        MUL(G,H,G)
        if bit == 1:
            APPLY(12*a,T,J)
            MUL(G,J,G)
        a := 2*a+bit                         # code generation
    MUL(G,Z2,ROOT)
    APPLY(1,ROOT,H)
    MUL(H,ROOT,CUBE)

The seven doublings and three added terms in the chain
each contribute one APPLY and one MUL. In total this
sequence has fourteen of each, including the final cube
check calculation. It costs 14*(6641+21445)+32.
Initialize a flag accumulator to zero. For each h load
CUBE+h and KEY+h, XOR, and OR into the accumulator,
at four operations per word. Emit eq with zero and a
branch to EXTEND if the equality fails. These are 66
operations including the branch. If equal, emit
SCORE(ROOT,KEY), APPLY(zeta,ROOT,ROTATED),
SCORE(ROTATED,KEY), APPLY(zeta2,ROOT,ROTATED), and
SCORE(ROTATED,KEY). This costs 3*293+2*6641. Thus a
parity VISIT costs at most 407463 for every key.

In both traversals, subset tree depth and NEXT/LAST
have the following invariant. The depth is the number
of flips from the pass's initial received half.
NEXT[depth] is the next potential child coordinate.
LAST at a nonroot depth records the entering edge.
Emit the following traversal twice, using the appropriate
VISIT body and literal depth cap in each. The cap is 112 in
the systematic pass and 109 in the parity pass. In the first copy END_PASS is the
second pass initialization, and in the second it is FINISH.
Label names are local to their copy. Branches not taken
fall through.

    VISIT:
        execute the pass's VISIT action
    EXTEND:
        flag := eq(depth,pass_cap)
        branch if flag to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        flag := eq(i,1000)
        branch if flag to BACK
        j := add(i,1)
        store(adr,j)
        in the systematic pass, emit the save and update above
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        FLIP(i)
        jump to VISIT
    BACK:
        flag := eq(depth,0)
        branch if flag to END_PASS
        adr := add(LAST,depth)
        i := load(adr)
        depth := sub(depth,1)
        in the systematic pass, emit the restore above
        FLIP(i)
        jump to EXTEND

FLIP expands literally at both sites:

    word := shr(i,6)
    off := and(i,63)
    mask := shl(1,off)
    adr := add(KEY,word)
    t := load(adr)
    t := xor(t,mask)
    store(adr,t)

A new node starts NEXT at one beyond the entering
coordinate. On return, the last flip is reversed.
Consequently both passes enumerate each subset of
size at most their respective caps once, along a strictly
increasing path of coordinate indices. All flips are in 0..999
and stack indices in 0..112. FINISH literally loads
BEST+h and stores to 32+h for h=0..15, then halts.
The best message has zero padding because the flips
and all field macros preserve the 1000-bit
representation. The two persistent registers depth
and BestD, up to 33 field registers including the
sixteen H and seventeen T words within MONOMIAL, and
at most twenty reusable arithmetic and control
temporaries fit the 64-register limit.

Let N_j=sum_(i=0..j) binom(1000,i). A pass with cap j has N_j
visits, N_j-1 successful extensions, N_j failed extension
attempts, N_j-1 nonroot BACK episodes and one root BACK.
The first pass's extensions cost 21+98+3*562+226=2031
and nonroot BACK costs 13+98=111. Failed extensions
cost at most six and the root BACK two. Including
293 per VISIT, they cost at most 2441*N_112.
In the parity pass traversal overhead is at most
40*N_109 by the corresponding original counts
21, six, 13 and two, and its VISIT bound is 407463.
Thus that pass costs at most 407503*N_109. These
counts include all branches and jumps.

Initialization and output cost 33+96+2*33+32=227,
plus 28086 to initialize the systematic square and
cube. Therefore for every input

    W=2441*N_112+407503*N_109+28313

bounds the operation count. The exact exponent check
is W^8 <= 2^4105, certifying L=4105/8. The declared
computation checks the error certificate and arithmetic
identities, including the incremental update, and uses
the fixed program and count proof for the full search.
