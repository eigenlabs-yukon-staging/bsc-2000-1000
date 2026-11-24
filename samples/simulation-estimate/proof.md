# LDPC code: block error about 0.051, claimed at most 0.06

The code is the null space of the parity-check matrix printed by gen.py (fixed seed 1), restricted to a 1000-dimensional subcode, decoded by belief propagation with 50 iterations.

We simulated 1,000,000 frames on BSC(1/10) and observed 51,034 block errors. A Clopper-Pearson 99.9% upper confidence limit for the block error rate is 0.0517, so eps <= 0.06 with a wide margin.


## Decoding complexity (claim 1)

The decoder runs exactly 50 iterations with no early stop, so its operation sequence does not depend on y. Messages are binary64 log-likelihood ratios, and tanh and atanh are 256-entry preloaded tables indexed by a quantized value. The graph has 6000 edges. Per iteration and edge: the variable update is at most 8 ops, and the check update (a table lookup, a product over the 5 other edges and the inverse lookup) is at most 24. With address arithmetic and loop control at most 8 more, that is 40 ops per edge, 240,000 per iteration and 12,000,000 for 50 iterations. Initialization (2000 LLRs, about 10 ops each) and the hard decision and output (about 20 ops per bit) add under 100,000. So W <= 12,100,000 <= 2^24 = 16,777,216.
