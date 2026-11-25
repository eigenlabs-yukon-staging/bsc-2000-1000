# Nineteen-operation decoding on the BSC

All indices start at zero, with low bits first. Encode a 1000-bit
message u by sending u followed by another identical copy of u.
The first copy proves injectivity, so the encoder has exactly 2^1000
codewords in 2000 channel uses.

The decoder is the following straight-line bsc-ops-v1 program.
The instruction templates are expanded before execution. Each address
and the zero operand in the second template is a program constant.

1. For each literal i from 0 through 2, emit R = load(i), followed
   by store(32+i, R).
2. For each literal i from 3 through 15, emit store(32+i, 0).
3. Halt.

There is one register, no preloaded table, and no memory access beyond
address 47. The additions in the template addresses are performed at
expansion time. There are no runtime branches, calls or loops. The
program uses 3 loads and 16 stores on every input, totaling 19
operations. The last output word is zero, including all pad bits.
The exact exponent comparison with the manifest's L=p/q is
19^q <= 2^p, checked by the declared experiment.

## Average block error

The first 192 decoded bits equal the first 192 bits of the received
word, and the remaining 808 decoded bits are zero. The first 192 bits
match the transmitted message exactly when their channel noise is zero.
The other 808 bits match exactly when the corresponding message bits
are all zero. These events are independent since the uniform message
is independent of the channel noise. Thus the exact success probability
is S = 9^192/(10^192 * 2^808). No constraint on any of the other received
coordinates is needed.

For the reduced manifest fraction B=b/d, the exact comparison
1-S <= b/d is equivalent to
(d-b)*10^192*2^808 <= d*9^192 since all denominators are positive.
The experiment computes this inequality for the declared b,d.
The fraction is chosen by rounding toward larger error.

## Optimality within the nineteen-operation budget

This additional observation concerns any explicit encoder of the
required dimensions. Every decoder writes the sixteen output words
under the cost rule. With at most nineteen operations on every path,
it can read at most three received words, adaptively or otherwise.
All other memory initially contains only constants or zeros.
Condition on a full random tape for any randomized decoder, independent
of the message and channel. Computation, including any input-independent
state and all addresses chosen, is then determined by previous replies
to reads. Writes to an original input address cannot give information
about the original word's unread bits. Reads of previously observed
words also give no new information.

Reveal each newly observed genuine channel bit in order, within
a word and across successive word reads, over a binary decision tree.
A final input word can contain padding constants, which require no
revelation. The deterministic decoder can be simulated from these
revelations. For analysis, pad any terminating path by revealing
previously unread channel coordinates selected in increasing order until
exactly 192 bits have been revealed. The output is fixed on each
padded leaf. There are at most 2^192 leaves. Along a specified path all
its revealed coordinates are fixed and distinct. Conditional on any
fixed message the probability of following that path is at most
(9/10)^192, since the channel is independent across those coordinates
and each prescribed reply has probability at most 9/10.

For a leaf to contribute to success, the message must equal the
single output on that leaf. Uniformity gives probability 2^-1000 for
that message. Summing the bound on its conditional leaf probability
over leaves gives success at most 2^192 * 2^-1000 * (9/10)^192 = S.
The upper bound is independent of any fixed random tape, so it holds
after averaging over that tape. The displayed decoder attains S.
