# Repetition code, ML decoding with random tie-breaking

Encoder f(u) = (u, u). The decoder g is maximum likelihood: for each i, if y_i = y_(1000+i) it outputs that bit; if the copies disagree, both values are equally likely and it outputs a fresh fair coin. Bit i is decoded correctly when neither copy flips (0.81), or when exactly one flips and the coin is right (0.18 x 1/2 = 0.09), with probability 0.9. Bits are independent, so eps(f, g) = 1 - (9/10)^1000 exactly, and 9^1000 * 10^46 >= 10^1000 (exact_check) gives eps <= 1 - 10^-46.

## Decoding complexity (claim 1)

Copy the first copy: for j = 0, ..., 14, load word j and store it to word 32 + j (30 ops); load word 15, mask it with 2^40 - 1 and store it to word 47 (3 ops). This is a straight-line program, so W = 33 on every received word, and 33^20 <= 2^101.
