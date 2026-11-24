# Parallel Hamming planes and aligned direct-address components

## Encoder and channel layout

Use zero-based bit indices, low bit first. The first 768 message
bits form three successive groups of 256. Within one group,
view the payload as four 64-bit words X0,X1,X2,D. Write
A=X0 XOR D, B=X1 XOR D, C=X2 XOR D. Send A,B,C,D followed
by A XOR B XOR D, A XOR C XOR D and B XOR C XOR D.
Across the three groups this occupies input words 0 through 20.
At each of the 64 bit positions independently, the four data
bits and the three parity bits make a length-seven codeword of
a systematic Hamming code. Thus 192 component codewords carry
768 message bits in 1344 channel positions.

Next there are fourteen 16-bit components occupying the next
224 message bits. Each is sent as a 16-bit systematic value a
followed by its parity P(a), occupying input words 21..27.
Its low fourteen bits l use the XOR of MAIN_COLUMNS
at their one positions, and its high two-bit tag t
XORs TAG_OFFSETS[t] into that parity. The XOR of the
three nonzero tag offsets is 20151, so the tag rule
is nonlinear. The last eight message bits have
the analogous 16-bit encoding using TAIL_COLUMNS. Send this last
component at transmitted bits 1984..1999, in the low 16 positions
of input word 31. Set intervening transmitted positions 1792..1983
to zero. Their noise is allowed by the channel. The total length
is 1344+448+192+16=2000 and the total payload is
768+224+8=1000. The four data coordinates of a Hamming lane determine
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
It places the results in the twelve output registers
across the three groups.

## Component decisions

Consider first a main received value with low fourteen
bits l_r, tag b and high sixteen-bit parity z_r. Let
s=z_r XOR P(l_r). For a hypothesized tag t, correct the
low and parity parts to the nearest word with that
tag. After removing the offset TAG_OFFSETS[t] from
the received parity, this is the systematic linear
component with fourteen data and sixteen parity bits.
Its error syndrome is s XOR TAG_OFFSETS[t].

The function leaders builds a shortest error for each
syndrome in that linear component. Construct the
syndrome graph whose edges XOR a literal data column
or a single parity bit. A sequence of edge indices
with given target has an XOR error giving the same
target and with weight at most its length. Conversely
each error of a certain weight directly supplies a
walk of that length. Thus the least path length is
the least error weight. Breadth-first search visits
each syndrome first at minimum distance and records
the XOR error along that path. It has that distance
as its weight. It fixes ties by the traversal order
shown in construction.py. All syndromes are reached
because the sixteen parity-bit edges generate them.
The tail component uses the same rule with eight
data and eight parity bits.

For the main case, start with the chosen error for
each hypothesized t and add the tag distance wt(b XOR t)
to its weight. The least of these four weights equals
the full minimum distance from the received block to
a transmitted main word. Denote it D(s,b). Select the
smallest t attaining it. The low data of the selected
word is l_r XOR the fourteen low bits of the chosen
syndrome leader for that t. This defines main_decode
for every key, including ties. Every component has a
unique selected message by these rules.

## Fixed memory and main schedules

Write H=2^32-2^16. At every address from 48 to H-1
preload the main decoder's sixteen-bit output for
that address in the low word bits, higher bits zero.
At address H+r for r in 0..2^16-1 preload the tail's
eight-bit decoded value for received word r shifted
left by 32. The two ranges use the available words
and are disjoint. Their values are specified by
table_word in construction.py. Addresses 0..31
hold the input words and 32..47 are zero until final
stores. If a main key lies in the prefix or the top
stripe, its load has the indicated actual contents
and may produce extra bits after packing. The proof
counts success only on safe main keys.

After the 72 operations of the three Hamming plane
groups, for each literal input word 21..27 load it
and get its two 32-bit keys by AND and SHR. Read the
table at each key. Relative main key k goes to
accumulator 12+floor(k/4) and lane k mod 4. Move
the result to the accumulator at lane zero, or shift
left by the lane times 16 and OR it in for each
other lane. When each key is safe the result occupies
only its intended lane. There are seven input loads,
28 extraction and lookup operations and twenty
shift and OR operations for ten nonzero lanes.
The counted main total is 55.

Load word 31, OR with H, load the tail entry from
that address and OR it into U15, at cost four. Word
31 has zero bits above the low sixteen received bits,
so the lookup gives the specified tail. Then AND
U15 with 2^40-1 to clear padding on every input,
even with exceptional main keys. Store all sixteen
accumulators in order to words 32..47. The generated
schedule contains no runtime loop, call or branch.
It uses sixteen output accumulators and fewer than
twenty other registers. Dynamic addresses are either
one unsigned 32-bit main key or the valid tail
address. The count is 72+55+4+1+16=148.
The exact bound is 148^100<=2^721.

## Exact local success and block error

For a fixed main received word of quotient pair (s,b)
the joint probability of seeing that word and of its
ordinary decision being correct is
9^(32-D(s,b))/(2^16*10^32). Indeed the message must
equal its uniquely chosen decision, and then the
channel noise is fixed by the received word. Its
weight is the chosen minimum D(s,b). This argument
also holds when the nearest transmitted codewords are
tied.

For fixed (s,b), every choice of l_r gives exactly one
received parity z_r=s XOR P(l_r). Thus there are 2^14
received words with this quotient. Summing the above
numerator over all quotient pairs and multiplying by
2^14 gives the ordinary success numerator over
2^16*10^32. Exclude main keys below 48 and keys in the
stripe H..2^32-1. In the prefix the received tag and
parity are zero, and low data equals r. The subtracted
numerator is

    sum_(r=0..47) 9^(32-D(P(r),0)).

In the stripe the received parity is 65535. At fixed
quotient pair there are as many such keys as low data
words l satisfying P(l)=s XOR 65535. Sum that exact
multiplicity times 9^(32-D(s,b)) over all quotient
pairs. This gives the stripe numerator to subtract.
The prefix and stripe are disjoint. The resulting
safe main success probability A is therefore exact
for the event that the main key is in the safe region
and the ordinary per-block decision succeeds.

For the eight-bit tail let T be the sum of 9^(16-w)
over all its minimum syndrome-leader weights w,
divided by 10^16. Translation by a tail codeword
leaves its syndrome fixed and gives exactly one
selected received error per syndrome. Thus T is
the exact tail success probability with its regular
lookup.

The 192 Hamming events, fourteen main safe-success
events and the tail's success event concern disjoint
message bits and disjoint channel uses. Uniform
messages and independent channel noise imply that
the probability of their intersection is h^192*A^14*T,
where h is the exact Hamming success above. Every
decoded payload word is correct on that event.
The declared computation checks

    1-h^192*A^14*T <= 12499999999999992169698947/12500000000000000000000000

from the fixed data and exact counts, and checks the
tag XOR witness. It replays the fixed schedule
on sample channel words with the byte-identical
reference opmodel and compares against separate
component calculations and word packing. The all-input
cost and sufficient success event follow from the
program expansion and probability proof above.
