# Aligned direct-address product

Bits inside words and blocks are numbered from zero, low first.
The encoder splits the message into 62 blocks of 16 bits and one block
of 8 bits. A k-bit block a is sent as a followed by P(a), where P is
the width-k linear map whose columns are the literal MAIN_COLUMNS
or TAIL_COLUMNS in construction.py. That code's parity_table() fixes
P(a) by XORing precisely the columns at the one positions of a.
The layout carries 62*16+8=1000 information bits in 62*32+16=2000
channel uses. Every encoder block is systematic and hence injective.

## Component decoding and fixed memory

For a received block r=(a,v), its syndrome is P(a) XOR v. This is
linear and its kernel consists of the component codewords. The
matrix columns are the selected parity columns followed by the unit
vectors. For each syndrome choose the full error vector recorded
by leaders(), a breadth-first search from zero with edges taken in
increasing coordinate order and first discovery winning ties. Each
edge XORs one syndrome column, and its tracked error vector XORs
the corresponding unit error. The graph is connected because of
the final k unit columns. Every error vector of weight w gives a
path of length w by using its support, and a path of length w gives
an error vector of weight at most w by XORing its coordinate
toggles. Therefore the shortest path length is the minimum error
weight for the syndrome. Breadth-first search visits vertices at
their shortest path length, processing its FIFO queue from zero,
so the stored error is a minimum-weight leader and its weight equals
that length. This also fixes the tie choices on the realized codes.

The ordinary component decoder removes the chosen error and returns
the low k bits. Its corrected word is in the kernel and hence is the
unique encoded word for that output. On any fixed sent codeword the
ordinary decoder succeeds precisely when the channel error equals
the chosen error for the received syndrome. Thus ordinary success
is exactly the sum of the probabilities of the chosen errors,
one per syndrome. In length n a chosen error of weight w has
probability 9^(n-w)/10^n.

The actual decoder has a fixed read-only word at every address
48 through 2^32-1. Let H=2^32-2^16. At each address r in [48,H)
store its ordinary 16-bit component output, replicated in each
of the four 16-bit lanes of one word. At address H+r for
0<=r<2^16, store the ordinary 8-bit tail output for r shifted
left by 32. This is exactly table_word() in construction.py.
The prefix 0..47 keeps its specified input and output roles.
The entire address space has exactly 2^32 words. The table
formulas and small parity and leader arrays specify every value
offline without dependence on any received word. In particular
replication is part of the stored constant, with no corresponding
runtime arithmetic.

## Complete fixed RAM program

The compile-time generator schedule() gives the tuple sequence
of the complete bsc-ops-v1 program. String operands name registers,
integer operands are literals. Each tuple applies its one named
primitive and writes the indicated result, except store writes an
output address and move is a free register copy. All loops and
conditions in the generator expand based on literal positions.
The RAM program itself is straight-line.

One load obtains each input word j for j=0..30. Its two received
32-bit blocks occur in its low and high halves. A literal mask
extracts the first key, and a right shift by 32 extracts the
second. For each key the program loads directly from that word
address and ANDs the value with the literal mask for its output
lane. The output component with index i uses lane i mod 4 in
accumulator U[floor(i/4)]. The first lane in each accumulator is
copied by a free move. Later lanes are ORed into it. The mask
guarantees disjoint lanes on all inputs. In particular keys below
48 still cause exactly one RAM load, from an input word or a
zero-initialized output word. No store occurs until the end, so
those values have their original contents. Keys from H upward
read a tail table value. For either exceptional group the mask
still defines a 16-bit output, though the value may differ from
ordinary decoding. A main key in [48,H) does return the ordinary
component output in its designated lane, because of replication.

The last input word has its higher 48 bits zero by the input
convention. The program loads that word, ORs it with H to form
the tail address, loads the pre-shifted value there, and ORs it
into U15. The main blocks have used bits 0..31 of U15. The tail
value uses bits 32..39 and has all other bits zero. Stores then
write all sixteen output accumulators to words 32..47, leaving
the required high 24 padding bits zero. There are 16 accumulators
and 3 other registers Y,R,D. All addresses and values are legal
for all inputs, including the exceptional keys.

The program has 95 loads, 93 ANDs, 48 ORs, 31 right shifts and
16 stores, totaling exactly 283 operations, with no runtime
branch, call, return or index arithmetic. The declared experiment
counts the generated sequence and verifies 283^200<=2^1629
with exact integers. It also compares the reference-counter replay
against a separate output calculation for representative words.
The sample preload contains every exactly evaluated read-only
address accessed by that sample and the test detects any missing
entry. The full fixed table is defined by the formulas above.
The samples include main keys in the reserved prefix and in the
tail table. Neither the sample run nor table generation supplies
an all-input cost assumption.

## Exact success lower bound for the actual decoder

The ordinary main component has the following chosen-leader
histogram (weights zero through six):
(1,32,496,4960,28624,31103,320). Write A for the sum of count times
9^(32-weight). Its ordinary success is A/10^32. We keep only
ordinary successes whose keys lie in [48,H), where the actual
decoder agrees. The two disjoint excluded key groups can be
accounted for exactly under a uniform 16-bit component input.

For each fixed received prefix key r<48 the ordinary decoder
selects exactly one supported message, with its noise equal to
the selected syndrome leader. Thus the probability of ordinary
success and receipt of that r is 9^(32-weight(e_r))/(2^16*10^32).
The verified prefix loss numerator Q is the sum of these 48
integer powers. Here the received parity half is zero, so the
syndrome is P(r).

For the stripe r>=H the parity half of the received block is all
ones. Given a fixed leader e, ordinary success holds for any sent
component with channel error e. For the stripe to occur we need
the transmitted parity P(a) to equal 65535 XOR the high half of e.
The program counts the exact multiplicity of every parity value
among the 2^16 messages. Multiplying that target multiplicity
by 9^(32-weight(e)) and summing over all leaders gives the stripe
loss numerator S, again over denominator 2^16*10^32. This
calculation holds for the exact submitted map regardless of its
rank.

Each main component is therefore ordinarily correct with a safe
key with probability
a = (2^16*A-Q-S)/(2^16*10^32).
The tail uses its ordinary decoder on every input. Its leader
histogram is (1,16,120,118,1) and its exact success is
b = sum_w count_w*9^(16-w)/10^16.
The 62 safe-success events and the tail success event depend
on disjoint component messages and channel flips, which are
mutually independent under the required uniform message.
On their intersection every output component is correct.
Hence whole success is at least a^62*b even though the program
can read another input word on an excluded main key.
The declared experiment computes the histogram, stripe and
prefix losses, and both component probabilities exactly, and
verifies
1-a^62*b <= 499999967504042273/500000000000000000.
