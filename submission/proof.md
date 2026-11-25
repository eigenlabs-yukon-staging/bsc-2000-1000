# Thirty-seven-operation decoder with a protected tail

All positions are zero based and all words are little-endian. Use a binary
matrix H with 32 columns represented by bytes, in this order:

- The integers 1 through 15, in increasing order.
- The integers 16, 32, ..., 240, in increasing order.
- The two integers 17 and 34.

A coordinate with a unit column 2^j is called a parity coordinate.
There are eight such coordinates. The remaining 24 coordinates are the data
coordinates, in increasing order. For a 32-bit word v, H(v) is the XOR of the
columns in its support. The parity columns are a basis. For each 24-bit
data word d, there is a unique word v with those data bits and
H(v)=248. Encode d as this word v. This defines an affine 32-bit code
with 2^24 words.

## Overall encoder

For a 1000-bit message u, make a 2000-bit transmitted vector as follows.

- Put message bits 0 through 959 in channel positions 0 through 959.
- Put the affine encoding of message bits 976 through 999 in channel
  positions 992 through 1023.
- Put message bits 960 through 975 in channel positions 1984 through 1999.
- Put zero in every other channel position.

These groups of payload bits partition the message, and each is directly
recoverable from the corresponding uncorrupted channel positions.
The resulting encoder is injective, has exactly 2^1000 codewords, and
uses exactly 2000 channel bits.

## Table and decoder program

For every s in 0 through 255, let e(s) be the least integer among the
minimum Hamming weight 32-bit words e with H(e)=s. Every s occurs,
because H has rank eight. Given a 32-bit received word v, let
D(v) extract the 24 data bits from v XOR e(H(v) XOR 248).
The corrected word has H value 248, so D(v) is a definite payload
and ties are completely fixed. As the channel flip probability is
1/10, minimum Hamming distance selects a maximum likelihood codeword.

The fixed read-only table at addresses 48 through 2^32-1 inclusive
has the word D(a) shifted left by 40 at address a. This occupies
2^32-48 words, disjoint from input and output. Entries are defined
by the formula even if they have value zero.

For each literal i from 0 through 14, in order, emit

    R = load(i)
    store(32+i, R)

with both addresses as constants. Then emit these instructions.

    R = load(15)
    R = shr(R, 32)
    V = load(R)
    V = shr(V, 24)
    T = load(31)
    V = or(V, T)
    store(47, V)

Halt. After the first shift R is a 32-bit integer, so even when it is
below 48, load(R) is an in-range memory access. Its value can be any
word without invalidating output padding. The next shift clips that
word to 40 bits. The 48 padding bits of input word 31 are zero for every
received word, so the OR also fits in 40 bits. If R>=48, the upper
24 bits of that output field are precisely D(R) and its lower sixteen
bits are precisely T. For R<48 the program still defines a decoder
with valid output padding, and the error analysis gives those inputs
no success credit.

This program uses at most three 64-bit registers and no scratch
memory. There are no runtime loops, branches, calls or implicit address
computations. Fifteen pairs followed by seven operations cost 37
operations on every received word. The exact exponent test for the
manifest L=p/q is 37^q <= 2^p.

## Exact success certificate

A syndrome consists of two nibbles. Each nonzero pure nibble is a
column, so any syndrome is reachable by choosing at most two distinct
coordinates, one for each nonzero nibble. There is one syndrome of
minimum distance zero, exactly 32 of minimum distance one because the
columns are distinct and nonzero, and the remaining 223 have minimum
distance two. For a uniform random data word and the 32-coordinate
BSC, the affine decoder D succeeds with probability

    (9^32 + 32*9^31 + 223*9^30) / 10^32.

Indeed affine syndrome H(v) XOR 248 is the syndrome of the noise.
On each of the 256 syndromes the decoder is correct precisely when
the noise equals its selected leader. This yields one success
contribution per syndrome, regardless of the transmitted data word.

The physical table cannot supply D on the 48 received values below
48. Those values have support in the first six bit positions. Their
syndromes H(v) are in 0 through 7 because those first six columns
are 1 through 6. Thus H(v) XOR 248 has high nibble 15 and low nibble
in 8 through 15. No column has this pair of nibbles, and the syndrome
is nonzero, so its least weight is two. For each such received v,
the joint probability that v is received and the transmitted data
is D(v) equals 2^-24 * 9^30/10^32. There is exactly one encoded word
for D(v), and it differs from v by the chosen weight-two leader.
Subtracting these 48 disjoint success events gives a lower bound on
what the actual table program succeeds on in the affine group.

The other 976 payload coordinates match exactly if their respective
channel positions have no flips. That event has probability
(9/10)^976 and is independent of the affine group's data and
noise. Even on an exceptional affine address, no credit is needed
for success. The whole program therefore succeeds with probability
at least

    S = [2^24*(9^32 + 32*9^31 + 223*9^30) - 48*9^30]
        * 9^976 / (2^24 * 10^1008).

For the reduced manifest fraction b/d, the claimed error upper bound
is 1-S <= b/d. The declared experiment reconstructs H, the exact
minimum-weight histogram through weight two, and all 48 exceptional
syndrome weights. It substitutes these counts into the displayed
success bound, checks the error inequality by cross multiplication,
and checks the exact complexity power inequality.
