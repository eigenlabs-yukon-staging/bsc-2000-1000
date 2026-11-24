# Extended Golay product with direct table decoding

All bit offsets are zero-based, low bit first. Split the message into 83
consecutive 12-bit blocks and a final four-bit block. For each 12-bit block
interpret its bits as coefficients of a polynomial a over F2, and form
c = a G using ordinary polynomial multiplication over F2, where
G = 1 + x^2 + x^4 + x^5 + x^6 + x^10 + x^11. Send the 23 coefficients of c
followed by their XOR as bit 23. Concatenate these 83 24-bit blocks, then
send the four-bit tail twice, consecutively. The lengths are 996+4
message bits and 1992+8 channel uses. Multiplication by the nonzero G is
injective, since F2[x] is an integral domain. The tail appears explicitly,
so there are exactly 2^1000 distinct encoded words.

## Chosen leaders and exact error

For a 24-bit word define its syndrome as the pair consisting of the
remainder of its first 23 bits modulo G and the parity of all 24 bits.
Both components are linear over F2. The kernel is exactly the image of
the component encoder. Indeed a divisible body of degree at most 22
has a quotient of degree at most 11, and overall parity zero determines
bit 23 uniquely. A received block therefore has the syndrome of its
noise, regardless of the transmitted component message.

For each syndrome choose the first error vector encountered when
enumerating weights zero through four, and for equal weights the
lexicographically increasing tuples of flipped positions. The declared
experiment computes the syndrome via polynomial long division, keeps
the first vector per syndrome and verifies this exact histogram by weight:

| Weight | 0 | 1 | 2 | 3 | 4 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Number | 1 | 24 | 276 | 2024 | 1771 |

There are exactly 4096 selected vectors, so every possible syndrome is
represented. No distance or weight-distribution theorem is needed.
The decoder XORs the chosen vector with the received block, obtaining a
kernel word, and recovers its unique 12-bit quotient. The tail output
copies only its first received copy. This determines every output.

For any transmitted block the decoder is correct exactly when the noise
equals the chosen vector of the corresponding syndrome. To see the
reverse implication from correctness, the corrected kernel word must
equal the encoded transmitted component message by injectivity. Taking
the XOR with the received word then recovers the actual noise uniquely.
Thus there is one successful noise pattern per syndrome. A specific
weight-w pattern has probability 9^(24-w)/10^24. If A is the integer sum
of histogram count times 9^(24-w), component success is A/10^24. Channel
coordinates for the 83 components and the four used tail positions are
disjoint, and each component's success criterion is independent of the
transmitted message. Consequently the conditional and average whole
message errors both equal 1 - A^83 * 9^4 / 10^1996. The experiment
computes A and this exact Fraction and compares it with the declared
B = 9999999874256283/10000000000000000.

## Fixed word RAM decoder

The direct table holds one word at address 48+r for each 24-bit
received value r. Its entry is the quotient defined above. The builder
fills the address for each pair of a component message and a selected
noise vector. This is well-defined and complete: the received syndrome
fixes the noise vector in such a pair, and removing it fixes the
component message. The number of distinct pairs is 4096^2=2^24. The
experiment materializes those entries and checks for duplicate writes.
Preloading is fixed before the decoder receives any input. The array
uses a compact host representation during validation, with one logical
RAM word per entry, so the address range including the input and output
prefix has length 48+2^24 < 2^32.

The complete program is the tuple sequence emitted by schedule() in
verify.py. Generation depends only on literal offsets. A tuple
(op,dst,args...) applies one listed RAM primitive to its operands and
assigns the result to dst. String operands denote registers. Store
takes its address and value from the operands and uses no destination
register. Move copies a register or constant at no cost. The loops
and conditions that generate the tuples execute at compile time. The
interpreter is a reference-counter driver for this straight-line
sequence. In the RAM program the only received-dependent memory
addresses are the table loads, and their additions and loads are
counted. Halting after the last store is free. Register use is 32 Y,
16 U and 3 temporary registers T, S and D, for 51 total.

The program loads the 32 input words into Y and initializes each U
register to zero. For a component beginning at off+64*j, shifting
Y[j] right by off brings the starting bit to position zero. If off
exceeds 40, fewer than 24 bits remain, and Y[j+1] shifted left by
64-off supplies the remainder with no overlap. At off=40 the shift
already leaves exactly 24 bits. All other cases mask to 24 bits.
The resulting integer indexes the table and produces D below 2^12.
For output, D is shifted left by the component's offset modulo 64 and
ORed into the current output word, or copied at offset zero. If its
starting offset exceeds 52, D shifted right by 64 minus that offset
starts the next word. This writes the low and high portions in exactly
their consecutive, disjoint positions. A direct assignment at offset
zero or to the next word replaces a zero register, because components
are emitted in increasing order without gaps or overlaps. The first
tail copy starts at offset 8 of Y31 and contributes four bits at offset
36 of U15. Final stores write words 32..47, leaving bits 40..63 of U15
zero. Hence the same decoder analyzed for error is implemented for
every received word, including those with many flips.

Counting the emitted primitive tuples gives 115 loads, 99 left shifts,
99 ORs, 83 adds, 83 right shifts, 74 ANDs and 16 stores, totaling 569.
Fixed scheduling gives this count on every input. The experiment checks
569^1000 <= 2^9153 with exact integers. Its sample runs with the
byte-identical reference counter independently exercise extraction,
table access and packing against the polynomial decoder. Those runs
illustrate the fixed program and do not replace the all-input argument.
