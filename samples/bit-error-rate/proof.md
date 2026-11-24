# Repetition code: error probability 1/10

Encoder f(u) = (u, u); decoder g(y) = (y_1, ..., y_1000), the first copy. Each decoded bit g(y)_i = u_i + Z_i is wrong exactly when Z_i = 1, which happens with probability 1/10. So the probability that a decoded bit is wrong is exactly 1/10 for every message and every position, and the error probability is at most 1/10.


## Decoding complexity (claim 1)

The decoder g (first copy), as a straight-line program in the bsc-ops-v1 model: for j = 0, ..., 14, load word j and store it to word 32 + j (30 ops, unrolled, no branches); load word 15, AND it with the free constant 2^40 - 1 and store it to word 47 (3 ops). Words 0..14 hold y_1..y_960, and the mask keeps y_961..y_1000 and clears bits 40..63, so the output is g(y). The program has no branches, so it executes 33 ops on every received word: W = 33. Since 33^20 <= 2^101 (exact integers: 33^20 is about 2.4 x 10^30 and 2^101 is about 2.5 x 10^30), W <= 2^(101/20).
