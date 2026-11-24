# Pruning shared parity tags with two coset lists

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

## Unrestricted comparison decoder

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

For a given tag t and received side (U,V), correcting with
l(q(U,V) xor J(t)) yields the nearest side word having this
tag. For comparison, define the unrestricted minimum-distance
decoder to minimize the sum of the corrections' weights over
all 1024 tags. At a tie it prefers any tag present in the two
lists defined below over any absent tag. It orders the present
tags by first occurrence during the actual list scan, and the
absent tags by numerical value. It uses the same private
leaders. The comparison tie rule therefore still gives a
minimum-distance decoder.

## Reduced tag lists

In the fixed table columns the low ten bits of J(t) equal t.
The declared runs check this identity, which makes J injective.
Define a twenty-bit quotient map of a thirty-bit syndrome x by

    Q(x) = (x xor J(x mod 1024)) div 1024.

Because J is linear and its low bits copy its input, Q is
linear with kernel exactly the range of J. For each h in the
twenty-bit quotient space choose, among the fifty-bit error
vectors e satisfying Q(q(e))=h, the first 64 in increasing
weight. Within a weight, use lexicographic order of the
increasing coordinate tuples of the flipped positions.
The all-tag side code is the kernel of Q composed with q.
Every quotient syndrome occurs since q is onto and Q maps
x=h*1024 to h. There are 2^30 errors per quotient syndrome
and the procedure defines each list completely.

Given a received side with q-value x, take its 64 chosen
errors e for h=Q(x). For each e include the candidate tag

    t = (x mod 1024) xor (q(e) mod 1024).

Then x xor q(e) belongs to the range of J since the quotient
syndromes agree, and t is the unique tag of the corrected
side. Repeated tags are allowed. The actual 100-bit block
decoder takes the union of the two sides' lists of tags,
scans side zero's list and then side one's list, each in its
construction order. It chooses the first tag with the least
sum of private-leader weights. It corrects
the two private words at the winning tag. This computes a
valid message on every received block.

## Exact error of the unrestricted comparison decoder

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
global occupied-volume bounds with these inputs give a
lower count D(t) of block cosets in the ball. The
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

The first quantity in the final bound is the probability lower
bound U0=N/10^100 for the unrestricted comparison decoder.
The restricted decoder agrees with a correct unrestricted
decision whenever the true tag occurs among the candidate tags.
The comparison decoder ranks the listed ties first, in the
same order as the actual decoder. If its winning tag occurs
in the lists, the restriction leaves the winner and the
private-leader decisions unchanged.

## Exact lower bound on side list hits

For each transmitted side word belonging to the all-tag
side code, its quotient syndrome is zero. Consequently
the syndrome of the received side equals that of the true
noise under Q composed with q. If this true noise is one
of its syndrome's first 64 errors, correcting by it would
include the true tag in that side's list. The probability
of this error-list event is identical for every transmitted
side word since it depends only on the noise. These two
error-list events in a block are independent because their
side noises are independent.

Let K=64. For each t let m(h,t) count errors of weight at
most t in the sphere for quotient syndrome h. The chosen
error lists contain exactly min(m(h,t),K) of these errors
for every h, since they choose by increasing weight.
For any integer k>=2K-1 and nonnegative integer m,

    min(m,K) >= K * (2k*m - m*(m-1)) / (k*(k+1)).

For m<=K the coefficient on the linear term is at most
one and the quadratic subtraction is nonnegative. For m>K
the inequality follows from (m-k)(m-k-1)>=0. Summation over
syndromes proves a lower bound on the cumulative list size
from the volume V of the 50-bit Hamming ball and its count
A of ordered distinct same-syndrome pairs. As before, A
is the sum of H(50,d,t) against the weight spectrum of
the all-tag side code over nonzero words. Choose
k=max(2K-1,1+floor(A/V)), round the rational lower bound
up to an integer, and take the maximum with the previous
lower count. Apply the same summation by parts argument
through t=22. This yields the exact error-list probability
lower bound S0=Nlist/10^50.

The actual component success probability is at least
U0 - (1-S0)^2 because a correct unrestricted result can
be lost only if neither side's error-list event occurs.
Those events are independent, as proved above.
The twenty original components use disjoint channel noises
and independent uniform input blocks. Their decoded outputs
also depend only on their own received bits. The full error
is therefore at most

    1 - (U0-(1-S0)^2)^20 <= 999790266569828209/1000000000000000000.

The supplied spectra.json records the zero-tag spectrum,
the four nonzero-tag square-convolution contributions,
and each cohort's contribution to the all-tag side
spectrum. A whole side spectrum is the sum of the 1024
local tag histograms because J is injective. Each sector
run enumerates the 2^20 private words at each of its
256 tag values, checks its capped-degree two-side
polynomial sum and its full all-tag side contribution.
The first sector also reproduces the zero-tag histogram.
Each run then checks the combined integer certificate for
both the unrestricted success bound and the side list hits.
Four sector successes establish every supplied array.

## Fixed tables and primitive program

Use A=2^30, Bc=2^31 and D=3*2^30. At A store the
weight words indexed by the thirty-bit private syndrome,
with value wt(l(r)). At Bc use the same sized bank
with the low twenty bits of l(r). At D store the 2^20
values A OR P(u) and at E=D+2^20 store 1024 values J(t).
Starting at 2^26 store an array of 2^20 rows of 64 words,
for the twenty-bit quotient syndromes. Its (h,j) word
is J(q(e) mod 1024) for the j-th chosen error e of h.
This read-only list array ends at 2^27 and the other four
arrays start at or above 2^30. All are disjoint and fit
in the 32-bit memory. The fixed generators outside model
runtime determine the entries completely.

The fixed iteration ranges below are unrolled at code
generation. Literal offsets and constants
are free, as are moves. For each block, extract its two
sequential 50-bit side fields from their received offsets.
Load the aligned word and shift right by the within-word
offset if nonzero. When the field crosses a word boundary,
also load the next word, shift left by 64 minus the offset
and OR. AND with 2^50-1. This uses at most six primitives
per side. Use an AND and right shift to get U and V.
OR U with D, load, and XOR with V. Call this register S,
which equals A OR q(U,V). These are five further primitives
for each side.

Move 61 to Best and zero to Winner. For each source side h=0,1, emit the
following seven primitives to set up its list.

    Low := and(S[h],1023)
    Address := add(E,Low)
    Shift := load(Address)
    S0p := xor(S[0],Shift)
    S1p := xor(S[1],Shift)
    Row := shr(Sp[h],10)
    Row := shl(Row,6)

Here Sp[h] is the literal S0p or S1p for the chosen source.
Shift equals J(Low). Row becomes 2^26+64*Q(q(U,V))
for that source's private syndrome. Then for each j=0..63
emit the following program body, with After[j] the literal
label of the next operation after the body.

    Address := add(Row,j)
    Offset := load(Address)
    Address0 := xor(S0p,Offset)
    Address1 := xor(S1p,Offset)
    W0 := load(Address0)
    W1 := load(Address1)
    Metric := add(W0,W1)
    Improve := ltu(Metric,Best)
    branch if Improve is zero to After[j]
    Best := Metric
    Winner := xor(Shift,Offset)
    After[j]:

By linearity Shift xor Offset equals the full J of the
candidate tag for the listed error. The two addresses stay
in the weight bank and Metric is the candidate distance.
The strict comparison keeps the first candidate at a tie.
There are nine counted primitives before every branch
target, including the branch. The free move and one XOR
occur only after a strict improvement. Best starts at 61,
all weights are nonnegative, and Best decreases by at
least one per improvement. The number of extra XORs
is therefore at most 61 per entire block. After both
scans Winner equals J of the chosen tag.

Use AND with 1023 on Winner for the chosen Tag, at a cost
of one. For each side XOR its S with Winner, XOR with
the constant A xor Bc to change table bank, load the
correction, and XOR it with U. These four operations per
side recover the exact private words.
Pack them and Tag into a 50-bit output block using shifts
20 and 40 and two ORs. Keep sixteen output registers
initialized to zero. If the block lies within a register at
its literal output offset, shift left and OR. If it crosses,
also shift the original block right by 64 minus the offset
and OR into the following register. At most four primitives
per block suffice. These disjoint intervals end at bit 999,
so padding is zero. Store all sixteen registers to words
32..47. Besides the output registers the program keeps
fewer than 24 scalars at once, within the 64 registers.

The complete all-input worst-case count is at most

    W = 20*(2*(6+5) + 2*(7+64*9) + 61 + 1 + 2*4 + 4 + 4) + 16
      = 25336.

The declared arithmetic checks W^1000 <= 2^14629 in
integers. The operation schedule counts the actual
two-list decoder used by the average-error proof.
