# A field-cube code with 130-flip systematic search

## Encoder and decoder

Let m=1000 and let K be F2[X]/(P), where
P=X^1000+X^53+X^45+X^34+1. Represent an element in the polynomial
basis, coefficient of X^i stored in bit i. The declared computation
checks irreducibility. Specifically, it iterates modular squaring
from X and verifies X^(2^1000)=X modulo P, and gcd one of P with
X^(2^500)+X and X^(2^200)+X. To see why this suffices, an irreducible
factor of degree d divides X^(2^s)+X exactly when d divides s.
In its quotient field of size 2^d, Frobenius has order d since
its d-th power fixes the field by the order of the multiplicative
group, and any smaller positive power cannot fix 2^d roots of
a polynomial X^(2^j)+X of lower degree. Fixing the residue of X
is equivalent to fixing that entire quotient. Thus all irreducible
factors of P have degrees dividing 1000. Every proper divisor
of 1000 divides 1000/2 or 1000/5. The two gcd checks exclude them,
so the degree-1000 polynomial is irreducible.

The encoder is f(u)=(u,u^3), with field cubing in K and the two
m-bit halves sent in order. The first half is systematic, so all
2^m messages have distinct encodings. On receiving y=(a,b),
the decoder considers exactly the messages v with wt(v+a)<=130.
It selects one minimizing wt(v+a)+wt(v^3+b), breaking ties
by the smallest integer v in polynomial-basis bit order.
The candidates are nonempty since a is one. The decoder is
deterministic. The operation proof below counts this exact
restricted search.

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

## Exact block error bound

Condition on an arbitrary channel error z=(z_I,z_P) and let
i=wt(z_I), j=wt(z_P), w=i+j. The transmitted message u is
still uniform. If i>130, use an error bound of one.
Otherwise u is a candidate. For any nonzero systematic
difference d the competing message is v=u+d. Its systematic
distance from the received first half is t=wt(d+z_I).
It is eligible exactly when t<=130. At t>w it cannot win.
For t<=w, let ell=w-t. A strict improvement occurs when
the parity difference lies in the ball about z_P of radius
ell-1. A tie in total distance occurs when it lies on the
shell of radius ell about z_P.

Of the uniform messages giving a tie for this d, exactly
half put the alternative first under the numerical tie rule.
Indeed u -> u+d pairs these messages without a fixed point.
The parity difference is unchanged under this swap, and the
fixed noise z and systematic distance t leave the tie result
unchanged. The smaller integer of the two is the true
message in exactly one of the pair. Thus for the union
bound this competitor contributes its strict-win
probability plus half its tie probability. This is the
average of the probabilities for the parity balls of
radii ell-1 and ell. In particular no independence between
competitors is assumed.

Write V=V_ell+V_(ell-1) and C=C_ell+C_(ell-1).
The hyperplane calculation, using uniform mass 2^(1-m)
per parity difference, bounds the average for a balanced
d by (V+C)/2^(m+1). For an unbalanced d it gives
2V/2^(m+1). At systematic distance t there are at most
N_t=binom(m,t) nonzero differences d, since d+z_I is a
bijection. At most min(N_t,Kbad) can be unbalanced.
As V-C is nonnegative, adding over them gives conditional
error bounded by the minimum of one and
     
    sum_(t=0..min(w,130))
       [N_t*(V+C)+min(N_t,Kbad)*(V-C)] / 2^(m+1),

with V,C evaluated at ell=w-t. This bounds an error of the
actual decoder, since whenever the true candidate is present
an erroneous selection must outrank it by this rule.

The declared experiment computes this bound for w<222
and uses one for w>=222. The two half error weights
are independent, each with mass
M_i/10^m where M_i=binom(m,i)*9^(m-i). It sums the
conditional bounds against their product masses, using
one at i>130. All binomial, coefficient, mass and
rational computations are exact integers or fractions.
The result is at most 7197791244648179/125000000000000000.
This is an average over the actual uniform message and
independent BSC noise for the fixed field code.

## Word representation and fixed data

A field element occupies sixteen consecutive words,
least significant word first, with only 40 live bits
in word 15. Use the input and output addresses from the
rules, and the following scratch ranges of 16 words each:
YS=48, YP=64, U=80, SQ=96, CUBE=112, BEST=128.
Use NEXT=144 and LAST=275 for two arrays of 131 words,
and MULT=406 for a 256-word multiplication scratch table
with sixteen rows of sixteen words. All are disjoint.
The fixed square table begins at SQUARE_TABLE=662 and
has 125 blocks of 4096 words. At address
SQUARE_TABLE+4096*j+16*d+h is word h of the field
value (d*X^(8*j))^2, where d is the byte interpreted as
a polynomial, 0<=j<125, 0<=d<256, 0<=h<16.
The definition is implemented by square_table().

There are two more fixed tables, REDUCE_ONE=512662 of
length two and REDUCE_FOUR=512664 of length sixteen.
For s equal to one or four, the entry at offset c
is the carryless product of c and
P_low=X^53+X^45+X^34+1, represented as one word, for
0<=c<2^s. Even at s=4 the degree is at most 56.
The highest used address plus one is 512680.
Preloaded data depends only on P. The only decoder
working memory is inside the scratch ranges above.

## Primitive field macros

In this and the next sections, constant-range loops with
compile-time indices emit a fixed sequence of primitives.
These are macros and never issue runtime calls. Assignment
from a literal or copying a register is free. The only
runtime jumps and branches are explicitly listed in the
search program. All masks, offsets and constant addresses
are compiled into the instructions.

SHIFT_s operates on sixteen field-element registers R0..R15
for fixed s equal to one or four. It multiplies the represented
element by X^s modulo P.

    c := shr(R15,40-s)
    for h from 15 down to 1:                   # emit in this order
        t := shr(R[h-1],64-s)
        R[h] := shl(R[h],s)
        R[h] := or(R[h],t)
    R0 := shl(R0,s)
    R15 := and(R15,2^40-1)
    adr := add(REDUCE_s,c)
    t := load(adr)
    R0 := xor(R0,t)

The descending order retains the original previous word for
the carry. The calculation removes the high bits above
degree 999 and replaces them using X^1000=P_low in the
field. Every reduction fits wholly in R0. The count is
1+45+1+1+3=51 primitives.

SQUARE moves zeros into the sixteen accumulators R.
For each literal byte index j from 0 through 124,
load from U+floor(j/8), shift right by 8*(j mod 8), AND
with 255, shift left by four, and add the literal base
SQUARE_TABLE+4096*j. These are five operations.
For each word h=0..15 add h to that dynamic address, load
the table word, and XOR it into R[h], costing three.
Finally store R into SQ with sixteen literal stores.
By linearity of squaring it computes the square of U.
The count is 125*(5+48)+16=6641.

MULTIPLY computes SQ times U and puts it in CUBE.
First load SQ into sixteen registers and store zero
and SQ to rows 0 and 1 of MULT, for 16+32 operations.
Apply SHIFT_1 three times in succession, storing the
result after each into rows 2,4,8 respectively. These
three stages cost 3*(51+16). For each of the other
eleven nonzero rows d, in increasing order, write the
XOR of the row for the least set bit of d and the row
for the remaining bits. There are two loads, one XOR,
and one store per word, with literal addresses. These
rows cost 11*64. This defines row d to be
SQ times the degree-at-most-three polynomial d.

Now move zero into all R registers. For each literal
nibble index j from 249 down to 0, load U+floor(j/16),
shift right by 4*(j mod 16), AND with 15, shift left
four and add MULT to get the row address. Apply
SHIFT_4 to R. At each of the 16 word positions add the
literal word index to the row address, load, and XOR
into R. By Horner's method this accumulates SQ times U.
After the last nibble store the sixteen words in CUBE.
The count is
16+32+3*(51+16)+11*64+250*(5+51+48)+16=26969.

## Distance test and exhaustive search schedule

Initialize YS and U with the received first 1000 bits.
For each literal j=0..14 load input word j and store it
to YS+j and U+j. For j=15 also AND it with 2^40-1.
This costs 49 operations. To initialize YP, for each
j=0..15 load words 15+j and 16+j, shift the first
right 40 and the second left 24, OR, and store to
YP+j, costing 96. The last word of YP already has its
higher 24 bits zero, since input word 31 has only 16
live bits. Initialize depth to zero and BestD to 2001
with free moves. Store zero in NEXT[0]. This is 146
total initialization operations.

At a candidate U, the following macro EVAL is inlined.
Apply SQUARE and MULTIPLY. Start distance Dist and
comparison flag Smaller at zero by free moves.
For each literal h=0..15, emit these instructions:

    a := load(U+h)
    b := load(YS+h)
    t := xor(a,b)
    t := popcount(t)
    Dist := add(Dist,t)
    a2 := load(CUBE+h)
    b2 := load(YP+h)
    t := xor(a2,b2)
    t := popcount(t)
    Dist := add(Dist,t)
    b := load(BEST+h)
    lt := ltu(a,b)
    same := eq(a,b)
    Smaller := select(same,Smaller,lt)

Scanning from low to high makes the most significant
differing word determine Smaller. This compares full
numerical messages. Use ltu(Dist,BestD), eq(Dist,BestD),
AND of that equality with Smaller, and OR with the
strict-improvement flag to get Take. Select the new
BestD according to Take. For each h load U+h and
BEST+h, select the former if Take else the latter,
and store the result at BEST+h. These are four
operations per h. Thus EVAL costs
6641+26969+16*14+5+16*4=33903.
BestD=2001 ensures the first candidate always replaces
the initial contents. Subsequent replacements implement
the stated tie rule.

The search traverses the subsets of flipped systematic
coordinates up to size 130, listing each in increasing
coordinate order along its path. At a node the current
depth is the number of flips, and NEXT[depth] is the next
coordinate eligible as a child that has not already
been visited. LAST at a nonroot depth records the
coordinate on its entering edge. The following labels
and instructions run immediately after initialization:

    VISIT:
        EVAL                                  # inline macro
    EXTEND:
        flag := eq(depth,130)
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
        FLIP(i)                               # inline macro
        jump to VISIT
    BACK:
        flag := eq(depth,0)
        branch if flag to FINISH
        adr := add(LAST,depth)
        i := load(adr)
        depth := sub(depth,1)
        FLIP(i)                               # inline macro
        jump to EXTEND
    FINISH:
        for h=0..15:                          # emit literally
            t := load(BEST+h)
            store(32+h,t)
        halt

Both FLIP occurrences are this sequence of seven primitives:

    word := shr(i,6)
    off := and(i,63)
    mask := shl(1,off)
    adr := add(U,word)
    t := load(adr)
    t := xor(t,mask)
    store(adr,t)

Any branch that is not taken falls through. At each
node its first NEXT is one greater than the last
coordinate (zero at the root). Before entering a
child the parent index advances, and a return restores
the parent's U by toggling the child's coordinate.
Thus every subset is enumerated once in its unique
increasing order, including the empty subset, and
no others are enumerated. Each FLIP gets 0<=i<1000
and cannot change any padding bits. The stack range
is depth 0..130, with NEXT values 0..1000. The two
words of persistent state depth and BestD stay in
registers. At most 32 other registers are needed for
the field-element arrays even if both macros use
separately named 16-word arrays. Branch, address,
and arithmetic temporaries can share at most twenty
additional named registers by reusing them between
macros. This is at most 54 registers, within 64.
The finish copies precisely the best message and its
zero padding from BEST and halts.

Let N=sum_(d=0..130) binom(1000,d). The traversal has
N visits and N-1 entering and returning edges. Each
successful extension costs 21, including both branch
tests, seven stack primitives, the flip and the jump.
There is one failed extension for each node at cost
at most six. Each nonroot BACK costs 13, counting its
test, branch, address, load, subtraction, flip and
jump. The root BACK costs two. This is at most 40N
overhead operations. The 146 initialization and
32 output operations add 178. For all y, therefore,

    W = 33943*N+178

is an integer upper bound on the operation count.
The experiment computes W and verifies W^200<=2^113597,
which certifies L=113597/200. This argument includes
every runtime branch and back-edge. The experiment
checks the equivalent window arithmetic on examples
but does not run the full search, and the count and
decoder correctness come from the explicit program
and its invariants.
