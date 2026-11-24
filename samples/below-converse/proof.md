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

The decoder g (first copy), as a straight-line program in the bsc-ops-v1 model: for j = 0, ..., 14, load word j and store it to word 32 + j (30 ops, unrolled, no branches); load word 15, AND it with the free constant 2^40 - 1 and store it to word 47 (3 ops). Words 0..14 hold y_1..y_960, and the mask keeps y_961..y_1000 and clears bits 40..63, so the output is g(y). The program has no branches, so it executes 33 ops on every received word: W = 33. Since 33^20 <= 2^101 (exact integers: 33^20 is about 2.4 x 10^30 and 2^101 is about 2.5 x 10^30), W <= 2^(101/20).
