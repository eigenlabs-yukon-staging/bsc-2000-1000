# Repetition code, first-copy decoding

Encoder f(u) = (u, u), each message bit sent twice: an injective map onto 2^1000 codewords of length 2000. Decoder g(y) = (y_1, ..., y_1000). g(y) = u exactly when Z_1 = ... = Z_1000 = 0, so eps(f, g) = 1 - (9/10)^1000 exactly, and 9^1000 * 10^46 >= 10^1000 (exact_check) gives eps <= 1 - 10^-46.

## Decoding complexity (claim 1)

    j <- 0
    loop: r <- load(j); store(32 + j, r); j <- j + 1; if j < 15 goto loop
    r <- load(15); r <- and(r, 2^40 - 1); store(47, r)

The loop copies words 0..14 to 32..46 with one load and one store each, 30 ops in all, and the tail costs 3 ops. So W = 33 on every received word.
