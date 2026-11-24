# Parallel Hamming planes and aligned direct-address components

## Encoder and channel layout

Use zero-based bit indices, low bit first. The first 512 message
bits form two successive groups of 256. Within one group,
view the payload as four 64-bit words X0,X1,X2,D. Write
A=X0 XOR D, B=X1 XOR D, C=X2 XOR D. Send A,B,C,D followed
by A XOR B XOR D, A XOR C XOR D and B XOR C XOR D.
Across the two groups this occupies input words 0 through 13.
At each of the 64 bit positions independently, the four data
bits and the three parity bits make a length-seven codeword of
a systematic Hamming code. Thus 128 component codewords carry
512 message bits in 896 channel positions.

Next there are thirty 16-bit components occupying the next
480 message bits. Each is sent as a 16-bit systematic value a
followed by its parity P(a), occupying input words 14..28.
The parity map XORs the literal MAIN_COLUMNS of construction.py
at the one positions of a. The last eight message bits have
the analogous 16-bit encoding using TAIL_COLUMNS. Send this last
component at transmitted bits 1984..1999, in the low 16 positions
of input word 31. Set intervening transmitted positions 1856..1983
to zero. Their noise is allowed by the channel. The total length
is 896+960+128+16=2000 and the total payload is
512+480+8=1000. The four data coordinates of a Hamming lane determine
its four payload bits by the stated invertible transform.
The other components are systematic, so the encoder has
exactly 2^1000 distinct words.

## Hamming-plane decoder

The parity-check columns in increasing coordinate order within
a length-seven component are the three-bit labels
(3,5,6,7,1,2,4). For any received component the syndrome is
the XOR of the columns at its one bits. It is zero on all
encoded words and equals the syndrome of the noise. Each
nonzero syndrome is a distinct column and selects its unit
error, while syndrome zero selects zero error. The decoder
XORs off that error and returns the four systematic bits.
The corrected component has syndrome zero, whose kernel is
exactly the set of encoded words because given the four data
bits each parity coordinate is uniquely determined. Thus for
any transmitted codeword correction succeeds precisely when
the actual error is the selected error. Selection succeeds on
the zero error and all seven single-bit errors. Its exact
component success is h=(9^7+7*9^6)/10^7.

The four payload words in a group are the linear functions
A XOR D, B XOR D, C XOR D, and D of the data planes.
For a single selected correction, the first of these
functions flips precisely for syndrome labels 3 or 7,
the second for 5 or 7, the third for 6 or 7, and the
last for 7. These conditions are exactly S0 AND S1,
S0 AND S2, S1 AND S2, and S0 AND S1 AND S2, respectively.
Thus there is no need to calculate the four individual
data-coordinate flip masks before transforming the payload.

The three 64-bit syndrome registers S0,S1,S2 carry the
corresponding syndrome bit at each position for all 64
components of a group. The schedule computes the three
linear payload functions, reuses A XOR D and B XOR D in
forming the syndrome, ANDs the required syndrome bit
planes, and XORs the resulting masks into the payload
functions. Scalar word logic acts independently on the
bits, so this implements the corrected payload on every
received lane. For each group the schedule emits seven
loads, three raw payload XORs, six additional syndrome
XORs, four ANDs and four final XORs, giving 24 operations.
It places the results in the eight output registers
across the two groups.

## Direct-address components

For either systematic component of width k (16 or 8), a
received block r=(a,v) has syndrome P(a) XOR v. Select the first
minimum-weight error for each syndrome by breadth-first
traversal from zero with each coordinate considered in fixed
increasing order. This is leaders() in construction.py.
The syndromes form a graph with edges XORing a single parity
check column. Unit parity columns make it connected.
A weight-w error gives a path of length w and conversely
a path of length w gives an error of weight at most w by
toggling its coordinates. Therefore graph distance equals
minimum weight. A queue explores increasing graph distance,
and the first discovery gives the fixed selected error.

Ordinary component decoding XORs off the selected error
and returns the systematic part. As the corrected word is in
the kernel and every kernel word is uniquely determined by
its data, ordinary success is exactly the event that noise
equals its selected full error. A fixed weight-w error in
n bits has mass 9^(n-w)/10^n.

Let H=2^32-2^16. At each address r>=48 and r<H, preload
the ordinary width-16 decoder output for received r, replicated
in each of four 16-bit lanes of a word. At address H+r for
0<=r<2^16, preload the ordinary tail output for received r
shifted left by 32. construction.py defines each
of these fixed values as table_word(). Direct table addresses
fit exactly the allowed range. Addresses 0..47 retain the
specified input and output roles. All fixed table data and
ordinary component rules are fully determined by the literal
matrices and the explicit leader search.

The actual RAM schedule computes the Hamming registers,
then for every input word 14..28 loads its two keys.
A mask and a right shift by 32 extract the low and high keys
respectively. For each key it loads the word at that address
and masks to its designated output lane. The first lane of
an accumulator is assigned by a free move and other lanes are
ORed in. The ordinary decoder output is obtained on keys
in [48,H). On a lower key the load instead addresses an input
or output word. All final stores occur later, so these areas
still have their original input or zero values. On a top-stripe
key the load uses a tail word. All such exceptional outputs
are defined by the same mask, and their possible failures
are handled in the bound. Word 31 holds the tail received
block with no nonzero padding. The program loads it, ORs
with H, loads that fixed word and ORs it into U15. It then
stores all 16 accumulators at their literal output addresses.
The tail fills exactly bits 32..39 of the last output word,
so output padding is zero for every input.

## Exact probability and worst-case cost

The ordinary main width-16 success is A/10^32, where A
sums 9^(32-weight) over its selected syndrome leaders. On
any ordinary success with a key in [48,H), the actual
decoder gives the same output. Compute the ordinary
success mass outside this interval as follows. For each
received value r<48 an ordinary successful pair uses exactly
the selected message and its corresponding error, so
under uniform data its contribution has numerator
9^(32-weight(e_r)) over 2^16*10^32. Its syndrome is P(r)
since its received parity is zero. Sum these 48 terms for
the prefix loss numerator Q.

For the top stripe the received parity half is 65535.
Fix one selected full error e and condition on that error.
Ordinary success holds regardless of the transmitted
component, and membership in the stripe requires
P(a)=65535 XOR the high 16 bits of e. Count its exact
multiplicity over the 2^16 systematic values and weight
by 9^(32-weight(e)). Summing gives the stripe loss
numerator S over the same denominator. The two regions
are disjoint. Therefore the probability of an ordinary
successful main component in the safe interval is
a=(2^16 A-Q-S)/(2^16*10^32). The tail always applies
its ordinary decoder, with success b given by its exact
weighted leader sum over denominator 10^16.

The 128 Hamming success events, thirty safe main
success events and tail success event depend on disjoint
message bits and channel errors. They are independent
since the message is uniform and channel bits are
flipped independently. All outputs are correct on the
intersection, irrespective of the dummy channel bits.
The declared computation computes h,a,b and checks
1-h^128*a^30*b <= 999999999999709414494893/1000000000000000000000000
in exact arithmetic.

The tuple generator is compile-time code. A string operand
names a register and an integer is a literal. Every tuple
except move executes one RAM primitive, with store having no
destination register. Every emitted instruction is executed
for every received word, with no runtime branches, jumps or
calls. The Hamming part costs 2*24=48, retaining the eight
output words in registers. The direct main part costs
15+15*6+22=127, from fifteen input loads, six operations per
two keys and twenty-two output ORs. The tail costs four and all
stores sixteen. Thus the exact count is 195. There are sixteen
accumulators and fewer than twenty reusable temporary registers.
The complete memory range is 2^32 words with every load address
in it, including prefix aliases. The declared reference-counter
replay checks the schedule against separately decoded bit
components, including prefix and top-stripe cases. It checks
195^125<=2^951 exactly. The count and correctness on all
inputs follow from the fixed schedule, not the choice of tests.
