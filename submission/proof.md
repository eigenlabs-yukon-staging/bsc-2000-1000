# Shared parity tags for pairs of short systematic codes

## Encoder

Number coordinates from the least significant bit. The first twenty
columns in construction.py define a linear map P from twenty bits to
thirty. The remaining ten define a linear map J from ten bits to
thirty. The declared computation checks all 1024 images of J are
distinct. Define the side encoding by

    E(u,t) = (u, P(u) xor J(t)),

where the private twenty-bit field u comes first and the thirty-bit
field follows. For each successive 50-bit message block let u0 and u1
be its first and second twenty-bit fields, and t its final ten bits.
Send E(u0,t) followed by E(u1,t). Since P and J are linear, the
100-bit block code C is linear. It is injective because the first
twenty bits of each side recover its private block and J is
injective. The twenty independent blocks occupy all 1000 message bits
and all 2000 channel bits.

## Decoder and its fixed tables

For a side received as (U,V), let the private syndrome of any side
word (u,v) be q(u,v)=v xor P(u). Thus the set C0 of private codewords
E(u,0) is the kernel of q. For a given syndrome s choose a minimum
Hamming weight vector of syndrome s. Construction.py fixes this
choice by breadth-first discovery from zero in the syndrome graph.
For each vertex the outgoing edges are considered in coordinate
order, first the twenty private columns, then the thirty parity
units. Every vertex is reachable using parity units alone. A shortest
walk to a vertex can be interpreted by flipping the walk's edges
from the zero error vector. It cannot repeat a coordinate, since
removing two uses of a column would keep the endpoint and shorten
the walk. Conversely an error of weight w yields a walk of length w.
The first-discovered record is therefore a minimum-weight error and
its weight is at most thirty. Write l(s) for this error.

For a candidate tag t the syndrome of the received side relative to
its translate E(u,t) is q(U,V) xor J(t). Add l of that syndrome to
the received side. The corrected side belongs to the translate, and
its private bits give the chosen u. Its distance to the received
side is the minimum over all u for that tag. The two sides are
independent under a fixed tag, so the decoder minimizes the sum of
the two weights over all 1024 possible tags. Choose the smallest tag
at a tie, and use the stated leaders within that tag. The full
decoder concatenates these decisions. Each decision is a
minimum-distance word of the corresponding 100-bit block code.

Use four read-only table intervals. The weight bank starts at
A=2^30 and the correction bank at Bc=2^31. Each has one word per
thirty-bit syndrome. Weight entries are wt(l(s)), and correction
entries are the low twenty bits of l(s). At D=3*2^30 place 2^20
words with entry u equal to A OR P(u). At E=D+2^20 place 1024
words with entry t equal to A xor Bc xor J(t). These four ranges
are disjoint and end below 2^32. All are above the input/output
prefix. The generator specifies table values before the decoder
runs. No table construction occurs at decoding time.

## Exact block error certificate

For any binary linear code of length n, the map to the cosets
partitions all 2^n possible errors into regions of equal size. A
minimum-distance decoder has the same minimal weight on the entire
received coset. If a coset has minimum weight w, each received
word in it contributes exactly one successful pair of message and
received word, of joint probability 9^(n-w)/(M*10^n) for M
uniform messages. A coset has M received words, so its aggregate
success contribution is 9^(n-w)/10^n. This calculation works for
the deterministic ties above. Let F(t) count cosets containing an
error of weight at most t. A nondecreasing lower bound on F(t)
through a cutoff gives a lower bound on the success numerator by
summation by parts, as follows. If D(t) is that lower bound and
D(-1)=0, the numerator is at least

    sum_(t=0..T) (D(t)-D(t-1)) * 9^(n-t).

Expanding this sum, each intermediate coefficient on D(t) is the
positive difference 9^(n-t)-9^(n-t-1), and the final coefficient
is 9^(n-T). Replacing D by F cannot decrease the expression.
The omitted true-success terms are nonnegative.

For any finite set of selected errors, let V be its size and
m(s) its multiplicity in a coset. If the number of ordered distinct
pairs within a common coset is at most A0, the number of occupied
cosets is at least the ceiling of V^2/(V+A0) by Cauchy-Schwarz.
For any integer k>=1 a second bound is the ceiling of
(2kV-A0)/(k(k+1)). Indeed for each nonnegative integer m, the
inequality

    k(k+1) * 1[m>0] >= 2k*m-m*(m-1)

is equality for k and k+1 on the positive branch and otherwise
follows from (m-k)(m-k-1)>=0. Summation gives the bound. Both
expressions increase with V for fixed A0 and fixed k. Thus they
also apply with a certified lower value v<=V. The computation
takes k=1+floor(A0/v), the maximum of both ceilings and the
previous lower count. All the volumes used contain zero and v>=1.

For a fixed nonzero word c of weight d in a length n code,
the number of ordered errors (e,e xor c) both at weight at most t
depends only on d. Choose the a positions of e on the support of
c. On that support the two weights are a and d-a. Off of it they
agree, so one can choose at most t-max(a,d-a) of n-d positions.
The exact number is

    H(n,d,t) = sum_a binom(d,a) *
               sum_(0<=b<=min(n-d,t-max(a,d-a))) binom(n-d,b),

with a limited to 0<=a<=d and max(a,d-a)<=t. It vanishes for
d>2t. Summing over all nonzero codewords by their weight counts
gives the pair upper bound, and is exact for the full ball.

First apply this formula in length 50 to the full local ball
and the spectrum of the private code C0. Apply the two
occupied-volume inequalities to obtain lower cumulative leader
counts D0(t). The actual number of local leaders of weight at
most t is precisely the number of private syndromes in the
ball because the chosen leaders minimize weight. A lower bound
on the number of pairs of local leaders with combined weight at
most t is given by convolving the increments of D0 with
themselves and taking the cumulative sum through t. To see this,
list actual leaders in increasing weight. The j-th true
weight is at most the j-th in the synthetic list with
D0(t)-D0(t-1) entries of weight t, for every j present in the
synthetic list. This gives an injection from synthetic ordered
pairs of admissible total weight to true ordered pairs.

For the 100-bit block, select just the errors whose two halves
are individually chosen private leaders and whose total weight
is at most t. The preceding convolution is a lower bound on the
size of this set. If two such errors share a global coset, their
difference is a codeword. If its tag is zero, each of the two
local differences has syndrome zero under q. Uniqueness of the
chosen private leader then forces the two halves to agree on
both sides, so the errors are identical. Every ordered
distinct pair thus has a nonzero tag difference. Count these
pairs from above by summing H(100,d,t) times the weight-d
spectrum of codewords having nonzero tag. This drops the two
local leader conditions when computing the upper count. The
global occupied-volume bounds with these inputs give the
claimed lower count D(t) of block cosets in the ball. The
cutoff is T=17 and all terms beyond it can be omitted.

For a given tag, let h_t(w) count fifty-bit side words E(u,t)
at weight w as u runs over the private words. The weight
polynomial of codewords with that fixed tag in a whole block
is the square of the polynomial with coefficients h_t,
because the two private words range independently. Thus the
nonzero-tag spectrum is the sum of these squares for the 1023
nonzero tags. Its coefficients are at most 2^50, so the
experiment uses base 2^51 with polynomial coefficients packed
as digits and squares by ordinary integer multiplication.
The declared calculation only needs degrees through 34 for
H. Truncation after degree 34 preserves those coefficients.

The supplied spectra.json records the zero-tag side spectrum
and four disjoint groups of contributions to the nonzero-tag
spectrum. The four declared runs enumerate the 2^20 private
words and each run's 256 tag images. In sector zero the run
reproduces the private spectrum explicitly. Each run reproduces
its own 35 coefficients of the nonzero-tag sum. Together the
four equalities establish all data used by the bound. Each run
checks the subsequent integer arithmetic for the combined data.
It obtains the success numerator N from the lower counts as
above. The twenty component messages are independent under
uniform input and use disjoint channel noises, so full
success is at least N^20/10^2000. The experiment verifies the
claimed rational B=999376364468155741/1000000000000000000 for one minus this lower bound.

## Primitive program and exact worst-case count

The following fixed iteration ranges are unrolled at code
generation. Literal data, including the 1024 fixed displacements
J(t), are determined by construction.py. Moves of register
values and constants are free. The only decoder steps are the
listed primitives. No branch or jump is required.

For each of the twenty blocks in increasing order do the
following. Extract each consecutive fifty-bit side using its
literal received offset. Load the aligned word and shift right
by the offset within that word if necessary. When the field
crosses the boundary, load the next word, shift left by 64
minus the offset, and OR. AND with 2^50-1. This costs at
most six operations per side. For each side separately, use
AND and right shift to obtain the private field U and parity
field V. OR U with D, load its table entry, then XOR the
result with V. Call this register S. Its high bits give A
and its low bits give q(U,V). These are five further operations
per side. Retain the two U and S values.

Move 61 into Best and zero into Tag. For t=0..1023 in
increasing order, unroll the following instructions.

    address0 := xor(S0,J(t))
    address1 := xor(S1,J(t))
    weight0 := load(address0)
    weight1 := load(address1)
    weight := add(weight0,weight1)
    improve := ltu(weight,Best)
    Best := select(improve,weight,Best)
    Tag := select(improve,t,Tag)

The index addresses only the weight bank and every total is
at most 60. The strict comparison implements the stated tag
tie rule. There are exactly eight operations for every t.

Add E to Tag and load that entry, costing two operations.
For each side XOR its S with the loaded word to give an
address in the correction bank. Load that entry and XOR it
with the corresponding U. That gives the decoded private
word in three operations per side. Form the packed fifty
output bits by shifting the second private word left twenty
and Tag left forty and combining the three by two ORs.
The two shifts and two ORs cost four operations.

Maintain sixteen output registers initialized to zero by free
moves. To insert the packed block beginning at a literal
bit offset, shift left by the offset within its word and OR
with that register. If it crosses the word end, also shift
right by 64 minus the offset of the original packed block
and OR that result with the next register. At most four
operations are used per block. The packed fields are disjoint,
and the only bits set in these registers lie at message bit
positions below 1000. Finally store all sixteen registers to
words 32..47, so the final padding stays zero even at ties
or on any other received word.

The input words can be loaded by their literal addresses and
the output registers coexist with fewer than 24 temporary
or persistent scalar registers. Thus there are fewer than
64 registers. The count for all received words is at most

    W = 20 * (2*(6+5) + 1024*8 + 2 + 2*3 + 4 + 4) + 16
      = 164616.

The verifier checks W^1000 <= 2^17329 in exact integers.
The count proves the bound without sampling the decoder.
