# Nonlinear tagged components and a 268-operation decoder

## Encoder

Number coordinates from zero, with the low bit first.
Split the 1000-bit payload into sixty-two consecutive
16-bit blocks and one final eight-bit block. A main
block a has low fourteen-bit part l and high two-bit
tag t. Write P(l) for the XOR of MAIN_COLUMNS at the
one positions of l. The block sends

    a | ((P(l) XOR TAG_OFFSETS[t]) << 16),

with the explicit constants in construction.py. The
tag offsets obey
32103 XOR 53530 XOR 58058 != 0, checked by the
declared experiment. In particular, for l=0 the two
single-bit tags 1 and 2 have parities whose XOR differs
from the parity of tag 3. This witnesses nonlinearity
of the transmitted component. The last block sends its eight payload bits followed by
the parity obtained by XOR of TAIL_COLUMNS over
its one positions. The layout occupies 62*32+16=2000
channel uses and is systematic, so all 2^1000 messages
are represented. The tags are treated as ordinary data
bits when determining uniform message probabilities.

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

## Fixed memory and scalar schedule

Let H=2^32-2^16. At every address r from 48 through H-1
preload main_decode(r) in the low sixteen bits of one
word, all higher bits zero. At address H+r with
0<=r<2^16 preload the tail decision at received
value r shifted left by 32. The rest of each of these
words is zero. This fully specifies each table word
via table_word. The two ranges are disjoint and use
the available addresses through 2^32-1. There is no
other working memory beyond the given input and output
prefix. In particular, addresses below 48 cannot
hold table data. Reads of 0..31 see the input words
and reads of 32..47 before final output see the
initial zero words.

The counted program is the compile-time fixed expansion
of schedule() in construction.py. For each literal input
word j from 0 through 30 emit a load. Extract its two
received main keys by a mask of its low 32 bits and a
right shift of 32, respectively. Use each key as a load
address. Map the key number k=2*j+high to output word
floor(k/4), with 16-bit lane k mod 4. When the lane is
zero, move the table word directly to that output
accumulator, with no charge for the move. Otherwise
shift left by 16 times the lane and OR into the
accumulator. The loop and conditions use literal indices
to choose instructions at code generation time.

For a key in 48..H-1 this inserts the correct component
decision in its own lane, with no effects on other lanes.
For a key outside that range its actual value is still
defined by the fixed memory and operations. A prefix
load may bring in other received bits. A main load from
the top stripe uses a word meant for the tail and may
set bits outside its own lane. The probability proof
below uses a sufficient event that excludes these main
keys. The operation bound and output format apply to
them as well.

After processing the 62 main keys, load input word 31,
OR with H, and use that address to load the tail value
from the top stripe. Input word 31 contains the final
16 received bits and zero high padding, so this obtains
exactly the tail lookup result. OR the result into U15.
AND U15 with 2^40-1 to clear all required high padding
on every input, including those with exceptional keys.
Finally store the sixteen accumulated output words
to addresses 32..47. Each accumulator has been assigned,
including U15 by main key 60. These stores are the
first writes to the output region.

There are 31 input loads for the main portion, two
primitives per main key for extraction and its
lookup, and two more for shift and OR for each of the
46 keys with nonzero lane. The four tail primitives,
one final padding mask and sixteen stores give
31+2*62+2*46+4+1+16=268. There are no runtime
branches, address computations apart from the listed
key operations, calls or loops. The use of sixteen
accumulators and three extra registers meets the
64-register limit. Every dynamic key is one unsigned
32-bit number, inside the memory bound on every input.
The exact exponent comparison is 268^1000<=2^8067.

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
lookup. The uniform input has independent payload
blocks, and all their channel noises are independent.
Consequently the probability of all sixty-two
safe main successes and the tail success is A^62*T.
On this event the actual scalar program packs the
correct payload lanes without interference, so its
block error is at most 1-A^62*T. The exact checker
verifies

    1-A^62*T <= 499999964602496721/500000000000000000.

It enumerates quotient distances, address multiplicities,
and syndrome leaders in host Python. It also replays
the fixed 268-operation schedule with the byte-identical
reference counter, checking padding and output against
a separate word packing calculation on sample inputs.
The complete worst-case count and success lower bound
follow from the expansion and probability argument above.
