# Baseline: repetition code, read the first copy

This baseline is deliberately weak. It shows the submission format and an exact computational step; it is not a competitive code.

## Code

- Encoder: f(u) = (u_1, ..., u_1000, u_1, ..., u_1000). Each message bit is sent twice, so f maps {0,1}^1000 to {0,1}^2000 and is injective. The code therefore carries exactly 2^1000 messages in 2000 channel uses.
- Decoder: g(y) = (y_1, ..., y_1000), the first copy. It is deterministic, outputs 1000 bits and does not depend on the noise.

## Error probability

Let y = f(u) + Z. Then g(y)_i = u_i + Z_i for i = 1, ..., 1000. So g(y) = u exactly when Z_1 = ... = Z_1000 = 0. The Z_i are independent Bernoulli(1/10), so for every message u

    Pr[g(f(u) + Z) = u] = (9/10)^1000.

Averaging over a uniform U gives eps(f, g) = 1 - (9/10)^1000 exactly.

## Bound

Claim: eps(f, g) <= 1 - 10^(-46), that is, B = (10^46 - 1) / 10^46.

The claim is equivalent to (9/10)^1000 >= 10^(-46), which is equivalent to the integer inequality 9^1000 * 10^46 >= 10^1000. The experiment `exact_check` (exact_check.py) evaluates both sides with Python's exact integers and prints the comparison. Its captured output is the evidence for this step. For orientation, log10((9/10)^1000) = 1000 * log10(0.9) ≈ -45.76, which is consistent with the exact check. This decimal value is not part of the proof.

Hence eps(f, g) = 1 - (9/10)^1000 <= 1 - 10^(-46) = B.

## Decoding complexity (claim 1)

The decoder g above, written as a program in the bsc-ops-v1 model of rules.md. The received word y is preloaded in words 0..31 (bit i of y is bit i mod 64 of word floor(i/64)); the output u goes to words 32..47 in the same layout, with bits 1000..1023 zero.

    for j = 0, ..., 14:   r <- load(j);  store(32 + j, r)        (15 x 2 ops; unrolled, no branches)
    r <- load(15)                                                  (1 op)
    r <- and(r, 2^40 - 1)                                          (1 op; the constant is free)
    store(47, r)                                                   (1 op)

The loop over j is fully unrolled, so the program is straight-line: it has no branches and executes the same 33 primitives on every received word. Total: 30 + 1 + 1 + 1 = 33 ops for every y, so the worst case is W = 33.

Correctness: words 0..14 hold y_1, ..., y_960 and are copied unchanged. Word 15 holds y_961, ..., y_1024 in bits 0..63; the mask keeps bits 0..39, which are y_961, ..., y_1000, and clears bits 40..63. So the output is exactly (y_1, ..., y_1000) = g(y), in the required layout.

Score: the claim is W <= 2^L with L = 101/20, which is equivalent to 33^20 <= 2^101. The experiment `ops_check` (ops_check.py) checks this with exact integers (and that 33^20 > 2^100, so L cannot be 5). It also runs the program on the counter `opmodel.py`, a byte-identical copy of the organizer's reference/opmodel.py (its sha256 is printed), on 22 received words: all report 33 ops and the correct output. Those runs illustrate the count; the worst-case bound is the straight-line argument above.
