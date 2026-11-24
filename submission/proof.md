# Extended cyclic code with two information-set scans

## Explicit code

There are 19 blocks that each encode 52 bits to 104 bits, and one
extended Golay block that encodes the remaining 12 bits to 24 bits.
We use low bit first in every block. Their concatenation has exactly
1000 information bits and 2000 transmitted bits.

For the wide block, let p=103. In F2[X], let G be the monic gcd of
X^p+1 with the polynomial whose exponents are the nonzero squares
modulo p. This determines G; the computation also gives its bit
representation 0xb1c29f41ef30b. The original 52 rows have bodies
X^i G for i=0 through 51 in positions 0 through 102, and parity bits
in position 103 making each extended row even. Denote their span by
C. Coordinates 0,...,102 are labeled by F_p and coordinate 103 by
infinity. Write I for positions 51,...,102 and J for the complementary
ordered positions 0,...,50,103. The declared experiment computes the
two systematic bases for C with respect to I and J by row elimination,
checking the identity restrictions in both cases. The encoder takes
its 52 input bits as the coefficients of the basis systematic on I.
This is injective. The tail multiplies its input polynomial by
G24=1+X^2+X^4+X^5+X^6+X^10+X^11, sends the 23 coefficients followed
by overall parity, and is injective as well. The entire code therefore
has 2^1000 distinct words.

## Distance certificate and spectrum of C

The experiment checks that the extended rows in the I-systematic basis
are pairwise orthogonal for the binary dot product and have weights
divisible by four. The independence furnished by the systematic
restriction makes C a dimension-52 self-orthogonal subspace in
104 coordinates. Thus C equals its dual. The weights of all its words
are divisible by four, since the support-intersection parity of any
two span elements is zero and wt(a+b)=wt(a)+wt(b)-2|supp(a) intersect
supp(b)|. The all-ones vector is in C because it is orthogonal to
every even-weight vector in C and C is self-dual.

Two coordinate permutations preserve C. They are x -> x+1, fixing
infinity, and x -> -1/x, exchanging zero and infinity. For each
permutation the experiment applies it to every systematic basis row
and checks span membership by exact reconstruction from the I
coordinates. Invertibility then gives equality of the permuted span.

These two permutations generate the determinant-one projective action
over F_103. The powers of the first give all matrices
U(a)=[[1,a],[0,1]], because the field is prime. The second matrix
S=[[0,-1],[1,0]] conjugates these to all lower transvections
L(a)=[[1,0],[a,1]]. Row additions and S eliminate the off-diagonal
entries of any determinant-one matrix, with S available for a zero
pivot. The remaining matrices diag(a,a^-1) for nonzero a are generated
as well, since U(a)L(-1/a)U(a)S=diag(a,a^-1). Consequently every
determinant-one projective permutation preserves C.

This action is transitive on unordered triples. An invertible
projective map sends any ordered triple of distinct points to
(infinity,0,1): send independent homogeneous representatives of
the first two points to the two coordinate axes, then rescale one
axis to obtain 1 for the third, whose coordinates are both nonzero.
If the map has nonsquare determinant, compose with x -> 1/x, which
preserves the target unordered triple and has determinant -1.
There is no element of order four in F_103^*, whose order is 102,
so -1 is nonsquare. The resulting determinant is a square, and
multiplying the matrix by a scalar makes its determinant one without
changing its projective action. Transitivity on triples implies
transitivity on unordered pairs and singletons by counting the
number of supersets of each size in an orbit of smaller sets. More
explicitly, an orbit of points or pairs would give a constant number
of its members in every triple. For points this forces all points to
have the same orbit by comparing triples with two shared points. For
pairs, comparing triples {a,b,x} and {a,b,y} gives that for any distinct
x,y, the sums of the indicator differences for ax versus ay and
bx versus by vanish for all distinct a,b outside {x,y}. There are
at least three such a, so each difference is zero. This forces the
pair indicators to be constant as well.

For any fixed support of weight w, take a uniform group element and
let X be the number of its moved points in I. For j=1,2,3, the expected
value of binom(X,j) is binom(w,j)binom(52,j)/binom(104,j), by
transitivity on j-subsets. These three factorial moments determine
the first three ordinary moments. They match the intersection
with I of a uniform w-subset of all 104 coordinates. That
distribution has mean w/2 and third centered moment zero by
exchanging I and J. Its variance is w(104-w)/(4*103), obtained by
summing the variances and pair covariances of w half-population
indicators sampled without replacement. The same mean, variance
and third centered moment therefore hold for X.

For weight four the variance is positive, so mean two implies that
some image has X at most one. For w=8,12,16 the variance exceeds one.
Set Y=X-w/2. Were Y always at least -1, the integer polynomial
Y^3-Y would be nonnegative and zero only at -1,0,1. Its expectation
is zero by the moment identities. The distribution would be
supported on {-1,0,1}, with variance at most one, a contradiction.
Thus each of these three weights has an image with X at most
w/2-2, hence at most six. A hypothetical nonzero word of weight
at most sixteen would have one of these four weights, and some
group image would have a nonzero systematic mask of weight at most
six. The experiment enumerates every such mask and checks that its
encoded weight exceeds sixteen. This proves minimum distance at
least twenty for C.

The experiment now derives the complete weight distribution by exact
linear algebra. If A_d is the number of C words of weight d, for
every w the character sum gives
2^52 A_w = sum_d A_d K_w(d), where
K_w(d) = sum_j (-1)^j binom(d,j) binom(104-d,w-j).
Indeed, for a binary vector z, sum_(c in C) (-1)^(c dot z)
equals 2^52 for z in C^perp and zero otherwise, by pairing with a
codeword with odd dot product in the latter case. Sum this identity
over z of weight w, then interchange sums and group the terms by
the size j of the overlap with a fixed weight-d word c. Self-duality
identifies the dual weight counts with A_w. Here A_0=A_104=1,
weights below twenty or above 84 other than the two endpoints are
zero (the latter by all-ones complementation), and weights not
divisible by four are zero. The exact system in the remaining
coefficients has a pivot in every unknown column and no
inconsistent rows, as verified by the experiment. Its unique
nonnegative integer solution is thus the distribution of this
particular code. No numerical code spectrum is imported.

## Wide decoder and probability bound

Define the admissible set E of vectors e of length 104 for which
min(wt(e_I),wt(e_J)) is at most six. For a received y, the decoder
finds the codeword whose difference e=y+c is in E and has least
total Hamming weight. It scans I first, then J, and within each
scan enumerates masks in increasing weight and then lexicographic
order of their increasing support tuples. A strict improvement
replaces the candidate; other candidates leave the earlier one
selected. For a set H equal to I or J, there is one codeword b
agreeing with y at H. Every candidate for the H scan is b plus the
codeword whose H restriction is a mask of weight at most six.
This exhausts precisely those differences with wt(e_H) at most
six, since H is an information set. At least one candidate exists
from the zero mask.

Let F_t be the number of distinct cosets of C meeting the subset
E_t of admissible vectors of total weight at most t. A coset here
is a translate of C in F2^104, with 2^52 members.
For each fixed y, the contribution of correct decoding to the
joint probability summed over messages is
2^-52 Pr[Z=y+c_selected]. The least admissible weight depends only
on the coset of y, because the possible differences for any y
are precisely its coset. All received words in that coset give
the same contribution even if their tied minimizers differ.
Summing over the 2^52 received words per coset shows that the
average success is exactly the sum of least-admissible-weight
noise probabilities, one term per coset. Equivalently it is
sum_t (F_t-F_(t-1)) 9^(104-t)/10^104, with F_(-1)=0.

The volume V_t of E_t is the sum of binom(52,i)binom(52,j)
over 0<=i,j<=52, i+j<=t and min(i,j)<=6. For distinct e,e' in E_t
in the same coset, their difference is a nonzero word of C.
For a given nonzero difference word of weight d, an upper bound
on the number of such ordered pairs is the size of the intersection
of the unrestricted radius-t Hamming balls about 0 and that word.
With i selected bits on its support and h off its support, both
weights are at most t exactly when h <= t-max(i,d-i). Multiplying
the corresponding binomial counts and summing gives the intersection
computed by full_intersection(). Multiplying by A_d and summing
d>0 bounds the total ordered distinct-pair collisions by M_t.
The case d>2t contributes zero.

Let m_s be the occupancy of a coset by E_t. Then sum_s m_s=V_t
and sum_s m_s^2 <= V_t+M_t. By Cauchy-Schwarz,
F_t >= ceil(V_t^2/(V_t+M_t)). Also every positive occupancy m obeys
1 >= m-binom(m,2), so F_t >= V_t-floor(M_t/2).
An integer refinement uses (m-d)(m-d-1) >= 0 for all integer m
and d>=1. For m>=1 this rearranges to
d(d+1) >= (2d+1)m-m^2. For m=0 the right side is zero, so summing
over occupied cosets and using the collision upper bound gives
F_t >= ceil((2d V_t-M_t)/(d(d+1))). Use d=1+floor(M_t/V_t).
The experiment takes the maximum of these bounds and the previous
bound at each t, starting from zero, through t=22. Denote the resulting
nondecreasing cumulative lower sequence by D_t. The value
A=sum_(t=0..22) (D_t-D_(t-1)) 9^(104-t) is a lower bound on the
numerator of average success. To see this, summation by parts
expresses this numerator bound as the sum of D_t times consecutive
decreases of 9^(104-t) for t<22, plus D_22*9^(104-22).
These coefficients are nonnegative, so replacing D_t by F_t can
only increase the result. Doing so counts all success contributions
through weight 22, and higher weights contribute nonnegatively.

## Tail and full block probability

For the tail, the syndrome consists of the remainder of the first
23 bits modulo G24 and the parity of all 24 bits. Its kernel is
the set of encoded tail words, since any divisible body of degree
at most 22 has a quotient of degree at most 11 and parity fixes
its extension. Choose one error per syndrome by increasing weight,
breaking ties by increasing lexicographic support. The experiment
enumerates weights through four and verifies the histogram
(1,24,276,2024,1771), accounting for all 4096 possible syndromes.
The decoder XORs the chosen vector into the received word and
outputs the quotient. It is right exactly when the actual error
equals the chosen vector, since each corrected codeword determines
its quotient uniquely. Consequently tail success is
T/10^24, where T=sum_w histogram[w]*9^(24-w).

The encoder and decoder separate all blocks. Uniform message
coordinates and independent flips make the component success
events independent under the required average. Whole success
is at least (A/10^104)^19*T/10^24. The declared experiment computes
these integers and checks the stated upper bound on one minus
that probability.

## Complete bsc-ops-v1 program

The following pseudocode is code generation for a word-RAM program.
All offsets, register names, and constants are fixed. Loops over
blocks, positions or output words in this pseudocode are unrolled.
The two explicitly indicated SCAN loops are the only runtime loops
in a wide component. All arithmetic in a primitive is word arithmetic
and shifts truncate to one word. A subscript on a register refers to
a literal register name. Moves and literal constants are free.

Memory begins with the required 48 input/output words. Next are two
delta tables in I,J order, with exactly 2N consecutive words each,
where N=sum_(i=0..6) binom(52,i)=23251684. For each mask in the
fixed ordering, the table holds the corresponding span word of
the relevant systematic basis, low 64-bit word first then its
high 40 bits. Finally a range of 2^24 words holds the tail direct
lookup table. Its entry at received r is the decoded quotient:
the preloader writes the quotient u at the address for
r=encode_tail(u) XOR e for each chosen e. This is complete and
unambiguous since the received syndrome identifies e, and its
removal determines u. All tables are fixed offline by the
explicit formulas and generators in construction.py. The total
memory prefix through these tables has 48+4N+2^24=109784000 words.
No dynamic memory is required.

FIELD(start,width) expands into these primitives, with start and
width literal and 1<=width<=64.

    j, s = divmod(start,64)                     # code generation
    R := load(j)
    if s != 0: R := shr(R,s)                   # code generation condition
    if s + width > 64:
        T := load(j+1)
        T := shl(T,64-s)
        R := or(R,T)
    if width < 64: R := and(R,2^width-1)

FIELD costs at most six primitives, and returns the requested
received bits in their low positions with all other bits zero.
The two word pieces have disjoint positions when joined. An OUTPUT
macro packs a k-bit component D (k=52 or 12) at literal output offset.
The 16 U accumulators are initialized to zero with free moves.

    j, s = divmod(output_start,64)             # code generation
    if s == 0:
        U[j] := D                              # free move
    else:
        T := shl(D,s)
        U[j] := or(U[j],T)
    if s + k > 64: U[j+1] := shr(D,64-s)

This costs at most three primitives. The two portions go into
disjoint output intervals for consecutive components; at word
boundaries any assignment without an OR initializes a word
unwritten by earlier components. This leaves the last output
word's high 24 bits zero.

For each wide block, FIELD extracts its low 64 and high 40 bits
into Y0 and Y1. Initialize BestW to 105, Best0 and Best1 to zero.
For H=I,J in turn, do the following. Move zero to B0,B1. For
each fixed coordinate in H, emit shr of its Y word by its local
offset, then and with 1 to form bit b. Use two select primitives
to choose the lower and upper word of the corresponding fixed
basis row, or zero, according to b. XOR those words into B0,B1.
These six primitives per coordinate reconstruct the word b_H
agreeing with Y at H. Set E0 := xor(Y0,B0), E1 := xor(Y1,B1).
Move the appropriate table base to Addr, and use the following
runtime loop, exactly 15 primitives per iteration.

    SCAN:
        D0 := load(Addr)
        Next := add(Addr,1)
        D1 := load(Next)
        N0 := xor(E0,D0)
        N1 := xor(E1,D1)
        C0 := popcount(N0)
        C1 := popcount(N1)
        Wt := add(C0,C1)
        Take := ltu(Wt,BestW)
        BestW := select(Take,Wt,BestW)
        Best0 := select(Take,N0,Best0)
        Best1 := select(Take,N1,Best1)
        Addr := add(Addr,2)
        More := ltu(Addr,End)
        branch if More to SCAN                 # else fall through

End is the literal first address after the table. The address starts
at the base and increases by two every iteration until End, so
exactly N iterations run, including the last not-taken branch. No
control path skips or repeats a candidate. The two scans update
the shared best registers and implement the claimed minimum with
first-encounter ties. After both scans emit two XORs of Y0,Y1
with Best0,Best1 to recover the selected codeword words, a right
shift of the low word by 51, a left shift of the high word by 13,
an OR, and an AND with 2^52-1. These six primitives recover
the I restriction in D. Expand OUTPUT with the component offset.

The tail uses FIELD on 24 bits, one addition for the tail table
base, one load of D and OUTPUT. At the end, emit one literal
store for each U accumulator to its required output address,
then halt. FIELD, the basis preparation, SCAN, and OUTPUT use at
most 25 distinct work registers even with separate names for
all shown temporaries and counts. Together with 16 accumulators
this is under the 64-register limit.

Each wide block costs at most
12+2*(52*6+2+15N)+6+3 = 649+30N.
The tail costs at most 6+2+3=11, and the final stores cost 16.
Thus W=19*(649+30N)+27=13253472238 bounds every received word.
There are no runtime calls, implicit loop tests or uncharged
address operations. The declared arithmetic verifies
W^500 <= 2^16813, giving L=16813/500. The worst-case bound is
established from the entire program; no sampled operation counts
are needed.
