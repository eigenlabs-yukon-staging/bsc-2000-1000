# Fixed field-cube product with a bounded systematic component search

## Encoder and decoder

Split the message into three successive blocks of 333, 333, and
334 bits. Each block represents an element of its field
K=F2[X]/(P) in the polynomial basis, coefficient of X^i in
bit i. Use the following two polynomials.
    
    P333 = X^333 + X^50 + X^47 + X^12 + 1
    P334 = X^334 + X^47 + X^9 + X^2 + 1

For each block of width m, encode u by (u,u^3) in its
field and concatenate the codewords in order. The systematic
halves give an injection from exactly 1000 bits into exactly
2000. Set R=47 in each component.

For each P, the declared calculation verifies X^(2^m)=X
modulo P and gcd one between P and X^(2^(m/p))+X for each prime
divisor p of m. For 333=9*37 the tested indices m/p are 111 and 9,
and for 334=2*167 they are 167 and 2. For any irreducible factor
of degree d, the residue of X generates its field, so it is
fixed by the j-th Frobenius power exactly when d divides j.
Frobenius to the d-th power fixes the field by the multiplicative
group order. A smaller positive power cannot fix the field by
the polynomial root bound. Therefore every irreducible factor
has degree dividing m, and any proper divisor of m divides
some m/p. The gcd tests exclude all proper divisors. A factor
of full degree m exhausts the monic polynomial P.

For each received component (a,b), consider every message v with
wt(v+a)<=R. Choose the message minimizing
wt(v+a)+wt(v^3+b), with the smallest numerical v winning ties.
The decoder concatenates these decisions. The set is nonempty
since it includes v=a. The primitive program below implements
this same decoder.

## Hyperplane collision bound for each component

Let Q=2^m. For a uniform u and a fixed d!=0, the parity difference
(u+d)^3+u^3 equals d^3((u/d)^2+(u/d)+1). In K, the map
z -> z^2+z is F2-linear with kernel {0,1}, and hence has image
of size Q/2. The trace Tr(z)=sum_(i=0..m-1) z^(2^i) is F2-linear
with values in F2 since it is fixed by squaring. It is nonzero
by the polynomial root bound, and Tr(z^2+z)=0. The parity difference is
uniform on the affine hyperplane Tr(d^-3 c)=Tr(1), where
Tr(1) is m mod 2. The relation between mask weights and the
absolute ball bias below applies to either right-hand side. Its nonzero coordinate mask
has entries Tr(d^-3 X^j). The map taking a field multiplier
to that vector of trace coordinates is linear and invertible:
a nonzero multiplier whose trace coordinates were all zero
would, by invertibility of field multiplication, make Tr
zero on all K. Each multiplier d^-3 has at most three
preimages by the degree-three root bound.

Let V_l be the m-bit Hamming-ball volume of radius l, and
let D=24. Masks of weight at most D or greater than m-D
are called unbalanced here. There are at most
Kbad=3*(V_D-1+V_(D-1)) values of nonzero d with these masks,
using complement for the high end. For a mask of weight s,
the signed sum (-1)^(mask dot e) over a centered radius-l
ball is, up to the sign of the center, the coefficient of
x^l in (1-x)^(s-1)(1+x)^(m-s). This follows by multiplying
the shell generating polynomial (1-x)^s(1+x)^(m-s) by
1/(1-x). For a balanced mask there are at least D pairs
of opposite factors. Multiplying each pair as 1-x^2,
then bounding coefficients in absolute value and using
1+x^2 <= (1+x)^2 coefficientwise on surplus pairs, yields
the bound C_l from (1+x^2)^D(1+x)^(m-1-2D). It satisfies
0<=C_l<=V_l. Accordingly the ball intersects the hyperplane
in at most (V_l+C_l)/2 entries for balanced masks and
V_l entries for unbalanced masks.

Let E_t consist of component errors of total weight <=t
and first-half weight <=R. Let J_t=|E_t|. For each possible
y let m_y count its candidates with errors in E_t.
Then sum_y m_y = Q J_t. For a fixed z in E_t, hold it
as the actual error and take u uniform. For a competing
u+d, its systematic error weight at y=f(u)+z is
s=wt(z_I+d), independent of u. It can count only for
s<=min(t,R), with allowed parity ball radius l=t-s.
The number of nonzero d at this systematic weight is
at most n_s=binom(m,s), and at most min(n_s,Kbad)
of them are unbalanced. Therefore the expected number
of competitors with errors in E_t is at most A_t/Q
for every fixed z, with the exact integers

    J_t = sum_(s=0..min(t,R)) n_s V_(t-s)
    A_t = sum_(s=0..min(t,R))
          [n_s(V_(t-s)+C_(t-s)) + min(n_s,Kbad)(V_(t-s)-C_(t-s))].

The second formula uses the mass 2/Q of each parity
difference in its hyperplane. Summing over every fixed
z in E_t and every true u bounds the count of ordered
distinct candidate pairs at every y by J_t A_t. Each
such pair at y occurs once as a first message and an
error z. Thus sum_y m_y^2<=J_t(Q+A_t).

## Coverage and block success

Let F_t be the number of y with m_y>0 divided by Q.
Cauchy-Schwarz gives F_t>=J_t Q/(Q+A_t). For integer
k>=1 use (x-k)(x-k-1)>=0 at integer occupancy x>0.
Rearranging yields k(k+1)>= (2k+1)x-x^2. At x=0
the right side is zero. Summing over all received
words gives the second bound
F_t>=J_t*(2k Q-A_t)/(Q*k*(k+1)).
Take k=1+floor(A_t/Q). The experiment floors both
bounds and takes the maximum with the preceding
lower bound, producing nondecreasing D_t<=F_t,
starting at D_(-1)=0.

For a received word whose least eligible distance is w,
the joint contribution of correct decoding summed over
uniform messages is Q^-1 9^(2m-w)/10^(2m). Exactly one
message is selected, and its noise probability is the
same for any tied minimizer. Summing by winning distance
gives exact component success
sum_t (F_t-F_(t-1))9^(2m-t)/10^(2m).
The experiment takes the partial numerator with D_t
through t=120. Summation by parts shows the resulting
A=sum_(t=0..120)(D_t-D_(t-1))*9^(2m-t) is no
greater than the actual success numerator. The
coefficients on D_t are positive differences of
successive weights and one positive final weight.
Later terms in the exact success are nonnegative.

Uniform component inputs and independent channel errors
make the component success events independent. Whole
success is at least A333^2*A334/10^2000, where each
A_m is the certified lower numerator for that width.
The declared computation checks the exact inequality
1-A333^2*A334/10^2000 <= 390125041029746637/1000000000000000000.

## Fixed data and RAM arithmetic

All offsets are low-bit first. Let k=ceil(m/64), h=m-64(k-1),
b=ceil(m/8) and n=ceil(m/4). In the field representation
there are k words and the high word has h live bits.
Initial memory uses the specified input and output prefix.
Use a separate bank for each field. The 333-bit bank starts
at 48 and is reused by the first two components. The
334-bit bank starts at the exclusive end of the first bank.
Starting at the bank base put six adjacent k-word arrays YS, YP,
KEY, SQ, CUBE, BEST in that order. The two arrays NEXT
and LAST follow, with R+1 words each. After them come
the 16 rows of k words each of multiplication scratch
MULT. The fixed reduction tables REDUCE_1 and REDUCE_4
have lengths 2 and 16 respectively and follow MULT.
At index c, REDUCE_s holds the carryless product
of c and P_low=P+X^m. Its degree is less than 64,
including at s=4, since deg(P_low)<=59.

After those tables put b consecutive square tables, one
for each byte position, with 256 k-word entries each.
The entry d in the j-th table consists of the field
words of (d*X^(8j))^2, with byte d treated as a polynomial.
These maps are fixed by the polynomial, and the supplied
table generators determine every word. The square tables are followed by SAVED with R+1 rows
of k words. AFFINE follows and holds m direction tables,
each with ceil(m/15) slices of 32768 entries of k words.
The field data module gives the same consecutive layout.
The exclusive bank ends are 1505885766 and 3016293468.
The full workspace and fixed table span fits within 2^32
words.

All fixed-range loops and literal conditions in the following
pseudocode expand at code generation time. A macro invocation
inlines its list, with no runtime call. Moves and literal
constants are free.

SHIFT_s, for literal s in {1,4}, multiplies the field value
in R0..R[k-1] by X^s. Emit the following sequence.

    c := shr(R[k-1],h-s)
    for j=k-1 down to 1:
        t := shr(R[j-1],64-s)
        R[j] := shl(R[j],s)
        R[j] := or(R[j],t)
    R0 := shl(R0,s)
    R[k-1] := and(R[k-1],2^h-1)
    adr := add(REDUCE_s,c)
    t := load(adr)
    R0 := xor(R0,t)

The descending order preserves the lower word's old bits
for each carry. The removed coefficients above m-1
are replaced by c*P_low in the low word, which performs
modular reduction. The length count is 3k+3.

SQUARE initializes k accumulators to zero. For j=0..b-1,
load KEY+floor(j/8), shift right by 8*(j mod 8), AND with
255, multiply the value by k, and add the literal base of
byte table j. This computes its dynamic entry address
in five primitives. For each word position add the literal
word index to that address, load, and XOR the value into
its accumulator. Store the accumulators to SQ at the end.
This costs b*(5+3k)+k. Field squaring is linear,
so SQ equals KEY squared. The high padding of KEY supplies
zeros if the final byte is partial.

MUL computes SQ times KEY. First load SQ into k registers,
writing zero and SQ to rows 0 and 1 of MULT. Use SHIFT_1
three times successively to write rows 2,4,8. For each
other nonzero row d in increasing order, XOR the stored
rows for the least set bit and the remaining bits and
store the result into row d, using two literal loads,
one XOR and one store for each word. This establishes row
d as SQ times the polynomial d. Move all R registers to
zero. For j=n-1 down to 0, load KEY+floor(j/16), shift
right by 4*(j mod 16), AND with 15, multiply by k and
add the literal base of MULT to compute the row address.
Use SHIFT_4 on R, then for each word add its literal index
to the row address, load and XOR into R. This implements
Horner's rule and ignores zero padding in a partial top
nibble. Store R to CUBE at the end. The count is
3k+3*(3k+3+k)+44k+n*(5+3k+3+3k)+k.
The non-runtime table preloads and the counted MULT
scratch updates are separate.

## Candidate program, padding, and operations

For each component, fixed code generation computes its
two input field starts. Extract every k-word representation
directly from the specified received words. For any
contiguous field piece of width w<=64 beginning at
offset s within word j, load word j and shift right
by s if nonzero. If s+w>64 load word j+1, shift it
left by 64-s and OR it into the current result. AND
with 2^w-1 when w<64. This is at most six primitives
and extracts precisely the desired low w bits with
zero padding. For each systematic piece use two literal
stores to copy it into YS and KEY. For each parity
piece store it into YP. This uses at most 15k
primitives. Move depth to zero and BestD to 2m+1,
and store zero at NEXT[0], costing one additional
primitive.

At initialization before the first VISIT emit SQUARE
and MUL once, writing the cube of the received systematic
field to CUBE. During traversal maintain KEY and its cube
using the edge operations below. The VISIT action zeroes
Dist and Smaller with free moves. At each literal word index j, emit

    a := load(KEY+j)
    b := load(YS+j)
    t := xor(a,b)
    t := popcount(t)
    Dist := add(Dist,t)
    a2 := load(CUBE+j)
    b2 := load(YP+j)
    t := xor(a2,b2)
    t := popcount(t)
    Dist := add(Dist,t)
    b := load(BEST+j)
    lt := ltu(a,b)
    same := eq(a,b)
    Smaller := select(same,Smaller,lt)

The most significant different word determines Smaller
because the sequence goes from low to high. Compute
ltu(Dist,BestD) and eq(Dist,BestD), AND the equality
with Smaller, OR the result with the strict-improvement
flag, and select BestD according to this choice.
At each word position load KEY and BEST, select KEY
if the choice is true or BEST otherwise, and store
it at BEST. These are 18k+5 scoring primitives.
The initial BestD lets the first candidate win,
and the update selects exactly the required minimum
with numerical ties.

The complete runtime control for this component is as
follows, reached initially at VISIT. All following
branches fall through if their condition is false.

    VISIT:
        emit the inline VISIT action above
    EXTEND:
        flag := eq(depth,R)
        branch if flag to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        flag := eq(i,m)
        branch if flag to BACK
        j := add(i,1)
        store(adr,j)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        SAVE_AND_UPDATE(i)
        FLIP(i)
        jump to VISIT
    BACK:
        flag := eq(depth,0)
        branch if flag to OUTPUT
        adr := add(LAST,depth)
        i := load(adr)
        RESTORE()
        depth := sub(depth,1)
        FLIP(i)
        jump to EXTEND

FLIP is inlined twice and does the following.

    word := shr(i,6)
    off := and(i,63)
    mask := shl(1,off)
    adr := add(KEY,word)
    t := load(adr)
    t := xor(t,mask)
    store(adr,t)

On descent before FLIP, SAVE_AND_UPDATE saves the k CUBE
words in SAVED+k*depth. Emit mul and add for the row base,
and one add, load and store per word, costing 2+3k. Let e=X^i
and b be the current KEY. The cube change is

    (b+e)^3+b^3 = e*b^2 + e^2*b + e^3.

The linear column at bit j is e*X^(2*j)+e^2*X^j.
The fixed table generator enumerates them with terms
advanced by X^2 and X. Add zero columns until the width
is a multiple of 15. The entry affine_entry(i,l,d)
XORs the fifteen columns selected by d in slice l, plus
e^3 if and only if l=0. Store each entry as k field words.
Its word h is at AFFINE+i*AFFINE_STRIDE+l*32768*k+d*k+h,
with AFFINE_STRIDE=ceil(m/15)*32768*k.

For the update, form the direction base with mul and add
and initialize k accumulators by loading CUBE. For each
literal slice l, load its low word of KEY, shift right
by its within-word offset, and AND with 32767. If it
crosses, load the next word, shift left and OR before
the mask. The crossing test at code generation is
(15*l mod 64)+15>64. The last slice fits in the allocated
high word and receives its unused zeros from the padded
KEY. Next multiply the key by k and add the literal
slice offset and direction base. For each word emit add,
load and XOR into its accumulator. After all slices
store the k accumulators. Thus a descent update costs

    U = 2+2k+ceil(m/15)*(6+3k)+3*crossings.

The resulting CUBE corresponds to the child reached by
the subsequent FLIP. On ascent RESTORE loads the saved
CUBE using mul and add for the row base and add, load
and store per word, at a cost of 2+3k. The following
depth decrement and FLIP restore the parent. Thus CUBE
is correct at every VISIT.

This visits the root and every increasing subset of
coordinate flips of cardinality at most R exactly
once. The parent advances its NEXT before the child
begins at the next greater coordinate. The child
records its entering flip in LAST and undoes it on
return. Thus depth is always the current cardinality,
and all flips stay in 0..m-1. The emitted control
checks and transitions account for these operations
on every path.

At OUTPUT, for each k-word piece at BEST, load it
and pack it into the specified output interval for
the component. If its starting offset in its output
word is zero, store it directly there. Otherwise
shift it left by that offset, load the old output
word, OR, and store. If the live portion crosses a
boundary, shift the original piece right by 64
minus the offset, load the next output word, OR,
and store. This uses at most nine primitives per
piece. In the direct-store case, the word has no
earlier bits since pieces are processed in
increasing order. In every other case the live
positions are disjoint, so OR accumulates the output
without altering any earlier piece. The high bits of
the final word remain zero. Continue by fallthrough
to the next component's emitted initialization or
halt after the last. The code-generation loop over
components adds no runtime jump.

For N=sum_(i=0..R) binom(m,i), each component has N visits.
There are N-1 successful extensions of 21 primitives,
N failed attempts of at most six, N-1 nonroot returns
of 13 and a root return of two. The traversal overhead
is at most 40N. Initialization and packing together
cost at most 24k+1 plus SQUARE and MUL for the initial cube.
The update, save and restore costs are charged once for
each visit in the upper bound. For each component use

    W_m = (U + 2*(2+3k) + (18k+5) + 40)*N
        + 24k+1 + SQUARE_cost + MUL_cost.

The integer sum 2*W_333+W_334 equals
W=13423322060905493287338313779110058800000132956489566120677370.
The exact comparison is W^1000 <= 2^203063.
The field macros use at most 2k field registers, together
with two persistent registers and at most twenty reused
temporary registers, within the 64-register budget.
The experiment checks the exact exponent and arithmetic
and exercises the multiplication and affine update identities on examples.
The all-input decoder count follows from the explicit
program, without running all its candidates.
