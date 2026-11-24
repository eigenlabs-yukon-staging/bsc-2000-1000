# Twelve-bit affine tags and packed correction entries

## Encoder

Number coordinates starting with the least significant bit. The
first nineteen columns in construction.py define a linear map P
from nineteen bits to thirty-one. The last twelve define a linear map J
from twelve bits to thirty-one, with 4096 distinct images as checked by
the declared experiment. A side with private word u and tag t
is the fifty-bit word E(u,t)=(u,P(u) xor J(t)).

Partition a twelve-bit vector into six consecutive two-bit
pairs. On each pair (low,high), let G act by (high,low xor high).
This is a linear map with G^2(s)=s xor G(s), and G(s)
differs from 0 and s for every nonzero s. The checker verifies
these facts for all inputs. For a pair of
twelve-bit parameters (a,s) define the four tags

    (a, a xor s, a xor G(s), a xor s xor G(s)).

For each successive 100-bit message component use the first
four nineteen-bit fields as its private words, the next twelve bits
as a, and the final twelve bits as s. Send E(u,t) for each of its
four sides in tag order. The first nineteen bits of each side
give its private word. Subtracting P of that word from the
parity field and using injectivity of J determines its tag.
The first two tags recover a and s. Thus the encoder is linear
and injective and the ten independent 200-bit components use
all 1000 message bits and exactly 2000 channel bits.

## Decoder and tables

For any side word (u,v), define the private syndrome
q(u,v)=v xor P(u). The 2^19 private codewords E(u,0)
form its kernel. Denote by l(r) a fixed minimum-weight
error of syndrome r. The generator chooses these by
breadth-first discovery from zero in the syndrome graph,
with outgoing columns in the order of the nineteen private
columns and thirty-one parity units. Each vertex is reachable
using parity units. A minimum walk does not use any
coordinate twice because dropping the two equal steps
would shorten the walk while preserving the endpoint.
The coordinates of a minimum walk therefore form a
minimum-weight error. Conversely any syndrome error of
weight w gives a walk of length w. In particular every
leader has weight at most thirty-one.

A side received as (U,V) under a hypothesized tag t is
decoded by correcting with l(q(U,V) xor J(t)). This gives a
word in the tag translate at the smallest distance achievable
with that tag. For each component the decoder considers all
4096^2 choices of (a,s), computes the sum of these four
side distances, and chooses the smallest sum. In a tie use
the first choice in the order (a div 8,s,a mod 8). The
returned private words are those from the corresponding four
side corrections. This is a minimum-distance decoder of the
actual 200-bit component linear code. The full decoder
concatenates the ten results.

Set A=2^31 and D=2^30. At A put one array indexed by
31-bit private syndrome. Its word for r equals wt(l(r))
OR the low nineteen bits of l(r) shifted left by eight.
At D place 2^19 consecutive entries A OR P(u).
At E=D+2^19 put the 4096 values J(t), then 4096 values G(s).
Their exclusive end is PROFILE_BASE. Starting there place
four arrays of 4096 scratch words for scalar profiles and
then 25 arrays of 512 words for packed profiles. Their
exclusive end is END=1074303488, less than A. All arrays
are above the protected input/output prefix. Both the table
bank at A and all the earlier ranges fit in 32-bit memory.
The scratch arrays are fully written before their use.

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
on the number of ordered four-tuples of local leaders with
combined weight at most t is given by the fourth convolution
power of the increments of D0 and the cumulative sum through t. To see this,
list actual leaders in increasing weight. The j-th true
weight is at most the j-th in the synthetic list with
D0(t)-D0(t-1) entries of weight t, for every j present in the
synthetic list. This gives an injection from synthetic ordered
four-tuples of admissible total weight to true ordered four-tuples.

For the 200-bit block, select the errors whose four sides
are individually chosen private leaders and whose total weight
is at most t. The previous convolution bounds the size of this
set from below. The difference of two errors sharing a global
coset is a codeword. If all the difference tags are zero, each
side pair has the same private syndrome and must agree, since
only one leader per syndrome was chosen. Every ordered distinct
pair therefore has a nonzero difference in the two tag parameters.
There are two classes of these differences. For s=0 and a
nonzero the four tags are all a. For s nonzero the four tags
are distinct and either none or one is zero.

For the differences with four nonzero tag slots, drop the
local leader conditions and count pairs using H(200,d,t)
against the full spectrum of codewords with this kind of
tag. For the differences with one zero tag slot, the two
errors in that slot must be the same leader. If that leader
has weight w, there are at most min(binom(50,w),2^31) choices
for it. The remaining three side words of the codeword
difference form a 150-bit word. For weight d of this word
the possible ordered errors in those sides are bounded by
H(150,d,t-w). Multiply by the number of choices for the
common leader and sum over w from 0 to t. Count the three-sided
words separately for each admissible tag tuple. The zero-tag
slot is already accounted for by the identical leader
condition, so this count has no multiplier of 2^19 from
the possible private words in the omitted slot.

For each local tag x write h_x(z) for the local weight
polynomial of E(u,x) over all 2^19 private words. Marking
the total three-side or four-side weights as powers of z,
the tagged spectrum polynomials are products of the
corresponding h_x. For s=0, a nonzero contributes h_a^4.
For s nonzero the sets {0,s,G(s),s xor G(s)} are additive
four-element subspaces. They partition nonzero directions
into triples. Each one is enumerated once by taking the
least nonzero s in it. For that subspace, the four-tag sets
are its additive cosets, with a chosen as the least
member of the coset. Each unordered coset yields twelve
distinct ordered parameter pairs, from the four initial
a choices and the three nonzero s choices, all with the
same polynomial product. The coset containing zero yields
a three-side product of h_s, h_G(s), h_(s xor G(s)),
with the same multiplicity of twelve. All other cosets
yield the four-side products. Together with the constant
nonzero tags these account for exactly 4096^2-1 pairs.

These polynomial sums are computed by packing integer
coefficients into digits of base 2^78. To bound digits let
M=max(1,max_x h_x(1/2)). The evaluation at 1/2 of any product
of up to four side polynomials is at most M^4. Each sum and
each intermediate polynomial in the computation consists
of at most 4096^2 such products, including integer
multiplicities. Since coefficients are nonnegative, its
coefficient at degree d<=50 is at most 2^d*4096^2*M^4.
The checker evaluates the bound at d=50 with exact
denominator 2^50 for each h_x(1/2) and proves it
is below 2^78. Thus ordinary integer products implement
all relevant convolutions without carries into or between
the retained slots.
The needed degree is at most 2t, and t is capped at 25,
so each product and sum can be reduced above degree 50.
The synthetic fourfold convolution of private leader
counts uses base 2^125 and coefficients at most 2^124.

The four supplied spectra0.json through spectra3.json
encode the 4096 tagged side histograms in consecutive
groups of 1024 tags. They store unsigned little-endian
three-byte coefficients in a Base64-wrapped zlib stream.
The four declared runs recompute all 51 coefficients for
the respective tag ranges [0,1248), [1248,2496),
[2496,3744) and [3744,4096) by enumerating all 2^19
private words, and compare each histogram to the record.
The fourth run also uses the full verified
arrays to compute the two tag spectra and the pair
limits, then the occupied-volume lower values D(t) and
the success numerator N. The components have independent
message bits and channel noises, making overall success
at least N^10/10^2000. The computation checks that
1-N^10/10^2000 <= 124410706854418567/125000000000000000.

## Scalar program and maximum operation count

Expand all fixed iteration ranges below at code generation.
Addresses indexed only by their unrolled indices and other
literals are literal constants in the program. Moves are free.
Instructions in the conditional block are executed only
through the stated branch. The straight-line blocks are
ordered by their unrolled indices.

For each component, extract its four 50-bit received fields
from the fixed input positions. Load the aligned word and
shift right by the within-word displacement when nonzero.
When the field crosses the word end also load the next word,
shift its contents left by 64 minus the displacement,
and OR it with the first piece. AND with 2^50-1. This uses
at most six primitives per side. Use AND and right shift
to give U and V. OR U with D, load the table entry and
XOR it with V to give S. These five additional primitives
make the high bits of S equal A and the low bits its
private syndrome q(U,V). Keep U and S for each side.

For each side i and tag t, XOR S with the literal J(t),
load the leader word, use AND with 255 to keep its
weight, and store it in the scalar profile at that tag.
This is 4*4096*4 operations. For each
side i, group g=0..511, and lane permutation p, gather
eight weights with tag indices 8g+(l xor p), l=0..7.
Emit a load, shift left by 8*l, and OR into an accumulator
for each lane. Store the accumulator in the packed profile
address for (i,p,g). This uses 25 operations per packed word.
For i=0 only p=0 is needed. For i=1,2,3 use all eight
p values. The pack build uses (1+3*8)*512*25 operations.

Move 125 into Best and zero into Slot. Move the constant
125*R into Limit, where R=0x0101010101010101. Let
H=0x8080808080808080. For each group g=0..511, load
packed profile word (0,0,g) into V0 once. Then for
each s=0..4095 the code has the three fixed tag offsets
d1=s, d2=G(s), d3=s xor G(s). For i=1,2,3 emit a
load of the literal packed profile word at
(i,di mod 8,g xor floor(di/8)). Call the words V1,V2,V3.
Lane l holds the side i distance for tag (8g+l) xor di.
Emit three additions to sum V0 through V3 in Sum.
Each lane weight is at most 124 since the four leaders
each weigh at most thirty-one. Thus no addition carries
to the next lane.

Emit the following five primitives after the three loads
and three additions for this group and s.

    Temp := or(Sum,H)
    Temp := sub(Temp,Limit)
    Temp := and(Temp,H)
    Skip := eq(Temp,H)
    branch if Skip past this candidate group to the next emitted operation

When computing the subtract, a lane's minuend is 128
plus its score, and its subtrahend is Best<=125. Each
difference is in [3,252], with no cross-lane borrow.
The high bit remains set exactly when that score is
at least Best. Accordingly the branch skips only when
no lane can improve the incumbent. Otherwise fall
through to the following instructions, unrolled in
increasing l=0..7.

    Candidate := shr(Sum,8*l)
    Candidate := and(Candidate,255)
    Improve := ltu(Candidate,Best)
    Best := select(Improve,Candidate,Best)
    Slot := select(Improve,4096*s+8*g+l,Slot)

Then emit Limit := mul(Best,R) and fall through to the
next body. This updates the replicated value of Best.
A non-skipped group decreases Best by at least one.
There are at most 125 such groups in all, since Best
starts at 125 and never becomes negative. Thus their
total extra cost is at most 125*41. All groups pay
11 primitives for the three loads, three additions and
five tests/branch above. The one V0 load per outer g
costs 512 more. Because the strict update keeps the
earliest minimizer, the traversal and tie rule agree.

At output, use AND and right shift on Slot to recover
a and s. Add s to the G table base and load G(s).
Build the three other tags using three XORs. These seven
operations recover the tags and both output parameters.
For each of the four sides, add E to its winning tag
and load the table entry. XOR it with that side's S to index
the leader bank. Load the leader word, shift it right by
eight, and XOR with U. This costs six operations
per side and gives the winning message's private word.

Keep sixteen output registers initialized to zero. Insert
the four private words and the two twelve-bit parameters
into their literal output positions. For each field
emit a shift left by the within-word offset and OR with
that register. If the field crosses, also emit a right
shift of the original field by 64 minus the within-word
offset and OR with the next register. There are at most
four primitives per field. Every field has its required
width, and the sixth field of the last component ends
at bit 999. Output padding remains zero. After all ten
components, store the sixteen words at addresses 32..47.
The sixteen output registers, eight saved U/S registers
and fewer than 24 other scalars fit in the 64-register
file. No other workspace is used.

Consequently the worst-case bound in counted operations is

    W = 10*(44 + 65536 + 320000 + 512 + 512*4096*11
            + 125*41 + 7 + 4*6 + 6*4) + 16
      = 234599456.

The exact comparison W^500 <= 2^13903 is checked by the
declared computation. The decoder operation count follows
from the emitted primitive program for every input.
