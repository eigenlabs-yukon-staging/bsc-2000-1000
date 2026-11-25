# Thirty-two operations using the padded input word

All bit positions are zero based, with low bits first. For a message
u of length 1000, transmit in the indicated positions:

| Channel bit positions | Contents |
| --- | --- |
| 0 through 959 | Message bits 0 through 959 |
| 960 through 983 | Message bits 976 through 999 |
| 984 through 1983 | Constant zeros |
| 1984 through 1999 | Message bits 960 through 975 |

The four rows partition the 2000 channel uses. Every message bit
occurs exactly once, so the encoder is an explicit injection of
all 2^1000 messages. Positions 2000 through 2047 in the RAM input
are padding according to rules.md, outside the channel. These 48
padding bits are zero for every possible 2000-bit received word.

## Decoder program and worst-case cost

Expand the following templates into a straight-line program.
All addresses are constants in the expanded program.

1. For each literal i from 0 through 14, in order, emit
   R = load(i), then store(32+i, R).
2. Emit R = load(31), then store(47, R).
3. Halt.

There is one live 64-bit register. The program accesses only words
0 through 14, word 31, and output words 32 through 47. It uses no
preloaded table, address arithmetic, branches, jumps or calls.
It executes sixteen loads and sixteen stores on every input, so
W=32 and W<=2^5 holds with equality. This is the full worst-case
count, including all the output stores.

Output bits 0 through 959 are the first 960 received bits. The lower
sixteen bits of output word 47 are the received bits at channel
positions 1984 through 1999. The high 48 bits of output word 47 are
zero since the source word has exactly those 48 padding bits zero.
In particular message bits 976 through 999 always decode to zero,
and the required output padding at positions 1000 through 1023
is valid even for a received word of all ones.

## Exact block-error bound

The decoded message is correct if and only if the 24 constant-output
message bits of the actual input message are zero and none of the
976 channel positions copied by the decoder has flipped.
The first event has probability 2^-24 for a uniformly distributed
message. It is independent of the other event, which has probability
(9/10)^976 by channel independence. Hence the average block error is

    1 - (9/10)^976 / 2^24.

For b/d equal to the fraction in the manifest, the desired upper bound
on block error follows from the exact integer comparison

    (d-b) * 10^976 * 2^24 <= d * 9^976.

The declared experiment computes this comparison directly. The manifest
fraction is rounded upward to meet the decimal-digit limit on scores.
