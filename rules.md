# Best [2000, 1000] code on the binary symmetric channel: error versus decoding cost

Construct a code that carries 1000 bits in 2000 uses of a binary symmetric channel with crossover probability 1/10. Prove two upper bounds for it: on its block error probability, and on the number of operations its decoder executes in the worst case. This is a Pareto challenge: a submission ranks by both numbers, and lower is better in each.

## Problem and model

The channel is memoryless: it flips each transmitted bit independently with probability exactly p = 1/10. A submission specifies:

- an encoder f: {0,1}^1000 -> {0,1}^2000, any function (linearity is not required); and
- a decoder g: {0,1}^2000 -> {0,1}^1000, either deterministic or randomized using randomness independent of the message and the noise.

The error probability is the average block error probability

    eps(f, g) = Pr[ g(f(U) + Z) != U ],

where U is uniform on {0,1}^1000, Z has 2000 independent Bernoulli(1/10) coordinates, U and Z are independent, and + is bitwise XOR. A bit-error rate does not count: one or more wrong message bits is a block error. A proved bound on the maximal error probability over messages also bounds eps and is accepted.

The code must be explicit. Both f and g must be completely determined by the submission, by one of:

- data files in the submission, for example a 1000 x 2000 generator matrix or a parity-check matrix;
- a construction stated in proof.md that fixes every parameter, for example a named algebraic code with its field, evaluation points, shortening and puncturing; or
- a deterministic Python 3 standard-library program in the submission that outputs the code, with every seed fixed.

The decoder may be "maximum-likelihood decoding of this code, ties broken by <a stated rule>" or any other exactly specified rule, but its cost is scored (below): the decoder must be given as a program in the decoding-ops model, and its worst-case operation count proved. Encoding cost is not scored.

## Decoding-ops model (bsc-ops-v1)

The second objective counts the primitive operations the decoder g executes, on a fixed abstract machine. `reference/opmodel.py` restates the model as an executable counter (sha256 `9b440361c3fede634b8021705f6afd62e964c1d996045c91a16228acee296903`); where the two differ, this text governs.

- **Machine.** A word RAM with 64-bit words and at most 64 registers. Memory has 2^32 words (addresses 0 to 2^32 - 1), initially zero. Moving values between registers and using constants written in the program are free.
- **Input.** The received word y is preloaded at no cost in words 0..31: bit i of y (counting from 0) is bit i mod 64 of word floor(i / 64), least significant bit first. Bits 2000..2047 are zero.
- **Output.** The decoder halts with the message in words 32..47 in the same layout, bits 1000..1023 zero. Each store is counted; halting and reading the output are free.
- **Tables.** Fixed read-only data (generator or parity-check matrices, lookup tables, schedules) may be preloaded at no cost from word 48 on. Tables are part of the decoder: the submission must determine them completely, and they count toward the 2^32-word memory. A decoder whose tables or working memory need more than 2^32 words is not in the model.
- **Cost: every primitive costs exactly 1 op.** The primitives are: load, store; integer add, sub, mul, unsigned div and mod, and, or, xor, not, shift left or right by any amount, popcount, count trailing zeros, equality, signed and unsigned less-than, select (conditional move); binary64 add, sub, mul, div, less-than, min, max, abs, and conversion to and from an integer word; one branch or jump (conditional or not, taken or not, including every loop test and back-edge, and each call or return); and rand, which returns 64 fresh uniform bits. Nothing else is a primitive: exp, log, tanh, square roots, vector or SIMD instructions and wider words must be built from these.
- **Sequential.** All operations are counted, on one processor. Parallel hardware gives no discount.
- **Worst case.** W is an upper bound on the number of primitives executed on every received word y in {0,1}^2000, and for a randomized decoder on every outcome of its random bits. An average over the noise, an expected iteration count or the count observed on sampled inputs is not a bound. Every loop, recursion and iteration count must be bounded by proof for every y; a decoder that may not halt has no bound.

The counted program must compute exactly the decoder g that claim 0 analyzes, on every received word. An iterative decoder with a fixed maximum number of iterations is a different decoder from the same decoder without the cap, and claim 0 must bound the error of the capped one.

## Scores

There are two scores, both lower-is-better, and a submission is compared with others by Pareto dominance.

### Block error, B (claim 0)

The score is the exact rational B in the claim "eps(f, g) <= B", given as a reduced fraction. Lower is better.

Two exact reference values (reproduce them with `python3 -I reference/bounds.py`, which uses integer and fraction arithmetic only):

- **Converse, about 0.052548118208.** Every code with these parameters has eps >= this value under every decoder. A decoding region of volume v can capture at most the probability of the v lowest-weight error patterns, which is a concave, increasing function of v. The regions are disjoint, so their average volume is at most 2^2000 / 2^1000. The protected minimum score is 525481183/10^10 (0.0525481183), just above the exact converse. The verifier rejects any claim below it before review.
- **Random coding, about 0.058506088975.** This is the random-coding union bound for i.i.d. uniform codebooks with ML decoding. It shows that some code reaches about 0.0585, but it names no code, so it does not count as a submission (see below). It marks the scale that explicit constructions are measured against.

The maximum score is 1. The best possible score lies between the converse and whatever the best code actually achieves, which is at most 0.0585.

Rounding and exceptional values:

- The submitted fraction is the score. It is never rounded, tightened or replaced by a weaker bound the reviewer finds in the proof; if the proof does not establish the claimed B, the submission is invalid.
- Yukon displays and ranks the nearest binary64 value of B. Distinct fractions can therefore display alike: the baseline's B = 1 - 10^-46 displays as 1.0. Where displayed values tie, the exact fractions decide.
- B = 1 holds for every code and earns no credit beyond showing the format. There is no score for an invalid submission, and no score for a partial result, because partial results are not accepted.

### Decoding ops, L (claim 1)

The score is the exact rational L in the claim "on every received word, g executes at most W <= 2^L operations in the bsc-ops-v1 model", given as a reduced fraction. The statement names the integer W, and the proof checks W <= 2^L exactly, as W^q <= 2^p for L = p/q. The score is on a log2 scale so that cheap decoders (tens of operations) and exhaustive maximum-likelihood decoding (more than 2^1000 operations for a general code) fit on one axis.

- The minimum is L = 4: any decoder stores at least the 16 output words. The maximum is 4096.
- The claimed L is the score. It is never tightened to the reviewer's own count and never rounded; if the proof does not establish W <= 2^L, the submission is invalid.
- Yukon displays and ranks the nearest binary64 value of L; where displayed values tie, the exact fractions decide.

The baseline's decoder copies the first copy of the message in 33 operations, so its claim is L = 101/20 (2^5.05 is about 33.1).

A person confirms official results.

## What does not count

- **Simulation.** Monte Carlo or simulation estimates, confidence intervals or empirical error curves. Simulation may appear as a sanity check, but it never establishes the bound.
- **Existence or ensemble bounds.** Random-coding, DT, RCU or Gallager-exponent bounds over an ensemble, and statements such as "a random linear code drawn from seed s satisfies the bound with high probability", unless the proof establishes the bound for the specific fixed code.
- **Different parameters or quantities.** Bounds for other lengths, dimensions or crossover probabilities, for bit-error rate, or for a channel model other than BSC(1/10).
- **Unproved inputs.** Bounds that rest on an assumed or estimated weight distribution, density-evolution thresholds, asymptotic results without explicit finite-length constants, conjectures, or numbers quoted from papers that are neither proved in the submission nor reproduced by its declared computations. Published theorems may be cited when every hypothesis is checked for this code.
- **Unchecked floating point.** Floating-point arithmetic without a rigorous error bound. Numerical steps must use exact integer or rational arithmetic, or interval arithmetic whose rounding direction the proof justifies.
- **Unrun computation.** Computational steps that the submission's declared experiments do not actually carry out within the run limits below. A computation too large to run there cannot support the claim.
- **Malformed codes.** Encoders or decoders that depend on the noise realization, codes that do not carry exactly 2^1000 messages in exactly 2000 channel uses, and decoders that do not output 1000 bits.
- **Average or sampled decoding cost.** An expected number of iterations, an average over the noise, or op counts from runs on sampled received words, offered as the bound W. Runs on the reference counter may illustrate a count; they prove the worst case only when the proof shows that the decoder's sequence of operations is the same on every received word.
- **Uncounted work.** Loop tests, branches, address arithmetic, loads and stores left out of the count; primitives not on the list (exp, log, tanh, SIMD, wide words) charged as one operation; tables or memory beyond 2^32 words, or tables the submission does not fully determine.
- **Two different decoders.** An error bound for one decoder (for example, maximum likelihood or an uncapped iterative decoder) paired with an op count for another (a cheaper approximation, or one with a cap the error analysis ignores).
- **Unscored lemmas.** Partial results without both bounds for a complete explicit code and decoder are not accepted.

## Submission and review

Edit only `submission/`. In `submission/submission.json`, declare exactly two claims, in this order. Claim 0 names the code and decoder and states eps(f, g) <= B, with `score` B. Claim 1 names the same code and decoder and states that g executes at most W <= 2^L operations of the bsc-ops-v1 model on every received word, with `score` L. Both scores are reduced fractions, and each claim lists its `evidence` (the files holding its proof). Put the complete proofs in `proof.md`, the decoder as a program in the model (pseudocode in proof.md, a Python program against a byte-identical copy of reference/opmodel.py, or both), and the code data or generator program in further files. Limits: at most 64 files, 262,144 bytes per file and 1,048,576 bytes in total.

Computations the proof relies on are declared as `experiments`: Python 3 standard-library programs that read their JSON input on stdin and print JSON. The verifier runs them in a sandbox with no network, up to 1024 MiB of memory, 2 CPUs and 120 seconds per run. A submission may declare up to 4 initial runs, using at most 480 seconds together. The reviewer may request one further run of up to 120 seconds. Printed results become evidence, and the reviewer reads every source file. A computation must check what the proof says it checks: for example, it should compute a bound exactly from the submitted code data, not print a constant.

One protected organizer tool is available, `bsc_exact` (optional, not required for a score). Declare it under `checks` in submission.json with input `{"kind": "tail" | "pairwise", "parameter": k, "bound_numerator": "...", "bound_denominator": "..."}`. It computes exactly, and compares with the rational bound given:

- `tail`: Pr[at least t of the 2000 bits flip] = Pr[Binomial(2000, 1/10) >= t], for 0 <= t <= 2001;
- `pairwise`: Pr[Binomial(w, 1/10) >= ceil(w/2)], for 1 <= w <= 2000. This bounds the chance that ML decoding prefers one fixed codeword at distance w over the sent one, counting ties as errors.

It returns `passed` (value <= bound), the exact value as numerator over a power of ten, and a 20-digit decimal rounded up. A pass establishes only that arithmetic fact. It says nothing about the submitted code: the proof must still establish, for example, the minimum distance or the number of codewords at each distance that the bound uses.

An LLM reviewer judges both claims together under the shared atomic-claims, evidence and external-evidence rules. Every step must be justified, and a single unjustified step makes the submission invalid. Structural validation is not mathematical acceptance. The review is an LLM assessment, not a formal certification, and a person confirms official results.
