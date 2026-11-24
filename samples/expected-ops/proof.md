# Repetition code, first-copy decoding

Encoder f(u) = (u, u), each message bit sent twice: an injective map onto 2^1000 codewords of length 2000. Decoder g(y) = (y_1, ..., y_1000). g(y) = u exactly when Z_1 = ... = Z_1000 = 0, so eps(f, g) = 1 - (9/10)^1000 exactly, and 9^1000 * 10^46 >= 10^1000 (exact_check) gives eps <= 1 - 10^-46.

## Decoding complexity (claim 1)

The program computes g(y) word by word. For each of the 16 output words it aligns the matching bits of the two copies (at most 11 ops: loads, shifts, ors, xor and the store of the first-copy word). Then it visits each position where the copies disagree: while d != 0, take i = ctz(d), record i in a statistics counter, and clear the bit with d = d and (d - 1). That is 5 ops per disagreement (ctz, add, sub, and, branch) plus one final branch per word. The output is the first copy, so the program computes g.

The copies of bit i disagree when exactly one of Z_i, Z_(1000+i) is 1, with probability 2 x 0.1 x 0.9 = 0.18, so there are 180 disagreements on average. The count is 16 x 11 + 16 + 180 x 5 = 1092 on average, so W <= 1097 <= 2^10.1.
