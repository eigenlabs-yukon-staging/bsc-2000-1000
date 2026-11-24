# Repetition code, first-copy decoding

Encoder f(u) = (u, u), each message bit sent twice: an injective map onto 2^1000 codewords of length 2000. Decoder g(y) = (y_1, ..., y_1000). g(y) = u exactly when Z_1 = ... = Z_1000 = 0, so eps(f, g) = 1 - (9/10)^1000 exactly, and 9^1000 * 10^46 >= 10^1000 (exact_check) gives eps <= 1 - 10^-46.

## Decoding complexity (claim 1)

Preload the table T with T[v] = v for every v in {0,1}^1000, 16 words per entry. The decoder loads the 16 words holding y_1, ..., y_1000 and masks the last one (17 ops). It forms the address of T[y_1..y_1000] by treating those 1000 bits as the index (15 shifts and ors), loads the 16-word entry (16 loads) and stores it to the output (16 stores). The program has no branches, so it executes at most 17 + 15 + 16 + 16 = 64 ops on every received word.
