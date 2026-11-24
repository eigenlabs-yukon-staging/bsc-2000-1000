# A nineteen-operation decoder for a binary repetition code

## Encoder

Number bits and machine words starting at zero. Within a
machine word, the lowest indexed bit is the least significant
one. Identify each binary message with its integer m in
0 through 2^1000-1 by that convention. The encoder sends
m followed by a second complete copy of m, meaning the
2000-bit integer m+2^1000*m. The encoder is injective
because its first 1000 transmitted bits recover m.
The BSC flips the 2000 transmitted coordinates independently
with probability 1/10.

## Decoder and operation count

Use a single 64-bit register T. For each literal i in
0 through 2 in increasing order, emit the following two
instructions, with i and 32+i substituted as constants:

    T = load(i)
    store(32+i, T)

Next emit one instruction store(32+i,0) per literal
i in 3 through 15. The program then halts. The twenty-four
padding bits are zero because the full last output word
is stored as zero. The program uses one register and
addresses below 48. It is fully unrolled and all addresses
and the zero constant are literals. There are no runtime
counters, branches, address calculations or preloaded
tables. The count is three loads and sixteen stores, so
W=19 for every received word and W^125<=2^531 by exact
integer comparison. The output consists of the first
192 received bits followed by 808 zero message bits.

## Error calculation

Success occurs exactly when the first 192 channel flips
are all zero and the last 808 message bits all equal zero.
The first condition suffices for those 192 output bits
because the first encoder copy is systematic. The second
condition is necessary and sufficient for the last 808
output bits because their decoded value is the constant zero.
Uniformity of the 1000-bit message makes all 808 selected
message bits independent and unbiased. Their values are
also independent of the channel noise. The probability
of the first condition is (9/10)^192 and that of the second
condition is 2^-808. The events impose no constraints on
the other transmitted positions. Thus the exact average
block error E is

    E = 1 - (9/10)^192 * 2^-808.

The rational calculation in verify.py reads the submitted
score, computes E and checks that B=ceil(10^255*E)/10^255.
In particular B is an upper bound and is strictly below one.
It verifies the power bound and exercises the operation
sequence on sample inputs. The operation and probability
arguments above apply to every received input and all
messages. The Python replay uses the reference opmodel.py and host
loops to simulate the fixed sequence, which includes no runtime
loops on the counted machine.
