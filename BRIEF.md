# BSC(1/10) [2000, 1000] code: error versus decoding cost — reconstruction brief

Written to the "Challenge reconstruction brief" template, on 2026-11-21, and revised the same day into a two-objective Pareto challenge (block error and worst-case decoding operations) at Sreeram Kannan's request. This challenge is not deployed. It is built for the shared verifier at `Layr-Labs/llm-verifier@da594aad5166918d4dd7d71630e1e649521e1953`, in this folder: `rules.md`, `verifier.json`, `checkers/bsc_exact.py`, `checker-tests.json`, `reference/bounds.py`, `reference/opmodel.py`, `submission/` and `samples/`.

## 1. Problem and target

Find error-correcting codes that carry 1000 bits in 2000 uses of a binary symmetric channel, and prove two upper bounds for each: its block error probability, and the number of operations its decoder executes in the worst case. The intended improvement is the Pareto front: lower proved error at a given decoding cost, or cheaper proved decoding at a given error.

- **Channel:** memoryless BSC with crossover probability exactly p = 1/10. Every one of the 2000 transmitted bits flips independently with probability 1/10.
- **Participant chooses:** an encoder f: {0,1}^1000 -> {0,1}^2000 (any function; linearity is not required) and a decoder g: {0,1}^2000 -> {0,1}^1000. The decoder may be deterministic, or randomized with randomness independent of the message and the noise.
- **Required theorems:** (claim 0) eps(f, g) <= B, and (claim 1) on every received word y in {0,1}^2000, g executes at most W <= 2^L primitives of the bsc-ops-v1 model (section 4). In claim 0, where eps(f, g) = Pr[g(f(U) + Z) != U]. U is uniform on {0,1}^1000, Z has 2000 i.i.d. Bernoulli(1/10) coordinates independent of U, and + is XOR. This is the average block error: one or more wrong message bits counts as an error. A bound on the maximal error over messages also bounds eps and is accepted.
- **Explicitness:** f and g must be completely determined by the submission, in one of three forms:
  - data files, such as a generator or parity-check matrix;
  - a construction in the proof that fixes every parameter;
  - a deterministic standard-library Python program, with fixed seeds, that outputs the code.

  "Maximum-likelihood decoding, ties broken by <stated rule>" is an exact decoder, but claim 1 must count a program that computes that same g on every input.
- **Organizer fixes:** n = 2000, k = 1000, p = 1/10, the uniform message distribution, the block-error criterion, and the bsc-ops-v1 cost model.
- **Excluded variants:**
  - other n, k or p;
  - bit-error rate;
  - other channels;
  - bounds over an ensemble or an existence argument instead of one fixed code;
  - average-case or sampled decoding cost, and encoding cost (not scored).

## 2. Acceptance and evidence

- **Standard:** rigorous proof under the verifier's default atomic-claims policy. One unjustified step makes the whole claim invalid; nothing is rescued by narrowing it.
- **Success probability:** not applicable. The claim is a deterministic inequality about an exactly defined probability.
- **Partial results:** not accepted. `allow_unscored` is false, so every submission must claim both B and L for one complete explicit code and decoder. Acceptance is atomic over both claims.
- **Evidence:**
  - The proof goes in `proof.md`.
  - Code data or a generator program goes in further files.
  - Every computational step must be a declared experiment the verifier actually runs: standard-library Python, JSON on stdin and stdout, no network.
  - Numerical steps must use exact integer or rational arithmetic, or interval arithmetic whose rounding direction is justified.
- **External results:** published theorems may be cited when every hypothesis is checked for this code. A number quoted from a paper (a weight distribution, a threshold) must be proved in the submission or recomputed by its experiments. The reviewer cannot fetch references.
- **Errors that invalidate an attractive result:**
  - simulation or confidence intervals offered as the bound;
  - an existence or random-coding bound for "some code";
  - an assumed or sampled weight distribution or minimum distance;
  - density-evolution thresholds, or asymptotic results without explicit finite-length constants;
  - floating point without error control;
  - a computation too large to run within the limits;
  - a code without exactly 2^1000 messages in 2000 uses;
  - an expected, average or sampled op count offered as the worst case;
  - loop tests, branches, loads or stores left uncounted, or non-primitives (exp, tanh, SIMD) charged 1;
  - tables or memory beyond 2^32 words, or tables the submission does not determine;
  - an error bound for one decoder paired with an op count for another.
- **Judge policy:** one authoritative reviewer (`gpt-6-astra`, high reasoning effort, 2400 s total review budget), with at most one follow-up batch. There is no committee or audit reviewer, and no reuse of prior judgments. Promotion on Yukon is manual: a person confirms official results.

## 3. Ranked score vector

| name | direction | exact formula and units | allowed range | source |
|---|---|---|---|---|
| `block_error_upper_bound` | minimize | B in the proved statement eps(f, g) <= B; a probability, as a reduced fraction | 525481183/10^10 <= B <= 1 | claim 0 |
| `decoding_ops_log2` | minimize | L in the proved statement "g executes at most W <= 2^L bsc-ops-v1 operations on every received word"; log2 of operations, as a reduced fraction | 4 <= L <= 4096 | claim 1 |

The track is a Pareto track (`mode: pareto` in benchmark.json): neither coordinate is a tie-break for the other, and a submission is dominated only by one at least as good in both.

- **Score from the claim:** each score is the submitted value exactly (B from claim 0, L from claim 1). It is derived only after the whole claim is accepted, and it is never tightened, rounded or replaced by a weaker bound the reviewer finds.
- **Allowed range of B:**
  - The minimum, 0.0525481183, sits just above the exact sphere-packing converse 0.052548118208…, which no code under any decoder can beat. The verifier refuses any claim below it before review.
  - The maximum is 1, and B = 1 holds for every code, so it earns nothing.
- **Allowed range of L:** at least 4, because any decoder makes at least the 16 output stores; at most 4096. The log2 scale puts a 33-op decoder (L about 5) and exhaustive maximum-likelihood search over 2^1000 codewords (L about 1008) on one axis; a raw op count of 2^1008 would exceed the verifier's 256-digit limit on score numerators. The statement also names the integer W, and the proof checks W^q <= 2^p exactly for L = p/q.
- **Reference scale, not a threshold:** the random-coding union bound for i.i.d. uniform codebooks with ML decoding is 0.058506088975…. So the best achievable code lies between 0.05255 and 0.05851, but that bound names no code and is not a submission. `reference/bounds.py` reproduces both values with integer and fraction arithmetic in under a second.
- **Rounding and exceptional values:**
  - Yukon displays and ranks the nearest binary64 values of B and L. Distinct fractions can display alike: the starter's B = 1 - 10^-46 displays as 1.0.
  - Where displayed values tie, the exact fractions decide.
  - An invalid submission has no score.

## 4. Complexity model and detailed results

One complexity resource: **worst-case decoding operations** of g in the fixed model bsc-ops-v1, stated as an upper bound W <= 2^L that holds on every one of the 2^2000 received words (and every outcome of the decoder's random bits). Encoding cost is not measured. The model is defined in rules.md and restated as an executable counter in `reference/opmodel.py` (sha256 `9b440361c3fede634b8021705f6afd62e964c1d996045c91a16228acee296903`); the text governs where they differ.

- **Machine:** a sequential word RAM, 64-bit words, at most 64 registers, 2^32 words of memory initially zero. Register moves and program constants are free.
- **Input and output:** y is preloaded free in words 0..31 (bit i in word floor(i/64), bit i mod 64, LSB first). The decoder writes u to words 32..47 in the same layout, bits 1000..1023 zero; the stores count.
- **Tables:** fixed read-only data (matrices, lookup tables) may be preloaded free from word 48 on, must be fully determined by the submission, and count toward the 2^32 words.
- **Unit cost:** every primitive costs exactly 1: load, store; integer add, sub, mul, unsigned div and mod, and, or, xor, not, shifts, popcount, count trailing zeros, equality, signed and unsigned less-than, select; binary64 add, sub, mul, div, less-than, min, max, abs, int/float conversion; one branch or jump (every loop test and back-edge, call and return); rand (64 uniform bits). Nothing else is primitive.
- **Quantifiers:** worst case over inputs and coins. No averaging over the noise, no expected iteration counts, no amortization over many received words, no parallel discount.
- **Retained constants:** all of them. W is an exact integer, and L = p/q with W^q <= 2^p.

**How the reviewer checks claim 1.** The reviewer works through five checks:
1. The counted program computes exactly the g of claim 0. A capped iterative decoder is a different g from the uncapped one, and claim 0 must analyze the capped one.
2. Every step is a listed primitive, and loop control, address arithmetic, loads and stores are all charged.
3. Every loop and recursion bound is proved for every received word.
4. Tables and memory fit in 2^32 words and are determined by the submission.
5. W <= 2^L is checked with exact integers.

Runs on the reference counter (a byte-identical copy in the submission) are supporting evidence: they give the exact count on the inputs run, and they establish the worst case only when the proof shows that the operation sequence is the same for every received word (a straight-line or fixed-iteration decoder).

## 5. Executable tools

**Protected checker `bsc_exact`** (`checkers/bsc_exact.py`): optional, organizer-supplied, not required for a score.

- **Purpose:** exact rational values of the two code-independent probabilities most BSC bounds are built from, compared with a rational bound.
- **Input:** `{"kind": "tail" | "pairwise", "parameter": k, "bound_numerator": "<decimal>", "bound_denominator": "<decimal>"}`.
  - `tail`: Pr[Binomial(2000, 1/10) >= t], the chance that at least t bits flip, for 0 <= t <= 2001.
  - `pairwise`: Pr[Binomial(w, 1/10) >= ceil(w/2)], for 1 <= w <= 2000. This bounds the chance that ML decoding prefers one fixed codeword at distance w, counting ties as errors, so it holds for any tie rule.
- **Output:** `passed` (exact value <= bound), the exact value as a numerator over 10^n, a 20-significant-digit decimal rounded up, and a message. Malformed bounds return `passed: false` with a message.
- **Algorithm:** a direct sum of C(n, j)·9^(n−j) over 10^n, using Python integers.
- **Data dependencies:** the protected constants n = 2000 and p = 1/10 in the source. It never sees candidate files, because the verifier's checker contract mounts only the organizer's bundle and caps inputs at 16 KiB.
- **What a pass establishes:** only the arithmetic fact value <= B. It establishes nothing about the submitted code: minimum distance, codeword counts and the decoder's behavior must still be proved.
- **Known-answer cases (`checker-tests.json`):**

  | input | expected |
  |---|---|
  | pairwise w=3 vs 7/250 | pass, value 28/10^3 |
  | pairwise w=3 vs 27/1000 | fail |
  | tail t=222 vs 1/20 | fail, value 0.056160877358706406131 |
  | tail t=0 vs 1 | pass, value 10^2000/10^2000 |
  | numerator "01" | fail (non-canonical) |

**Reference counter `reference/opmodel.py`:** organizer-supplied, not a registered checker. A registered checker cannot see candidate files, and a single run cannot establish a worst case, so this is a library: a `Machine` class with one method per primitive that counts each call, enforces the 2^32-word address space and reads out the decoded message. Participants copy it byte-for-byte into `submission/` and run their decoder on it in a declared experiment; the printed sha256 lets the reviewer confirm the copy. The starter's run gives 33 ops on all 22 tested received words.

**Participant experiments:** untrusted supporting evidence that the reviewer reads in full. The limits are:
- at most 4 initial runs, 120 s each and 480 s together;
- one follow-up run of 120 s;
- 1024 MiB of memory and 2 CPUs, no network.

A full 1000 x 2000 generator matrix (250 KB of bits) fits in the submission limits: 64 files, 256 KiB per file and 1 MiB in total. An experiment may check, for example, that the matrix has full rank, which shows the code carries exactly 2^1000 messages. A small or partial computation never certifies a full-length bound.

## 6. Starter and distinguishing cases

**Starter (`submission/`), expected to be accepted with scores B = (10^46 − 1)/10^46 and L = 101/20:**
- Repetition code f(u) = (u, u), decoded from the first copy, so eps = 1 − (9/10)^1000 exactly.
- The claim is B = 1 − 10^-46. The experiment `exact_check.py` evaluates 9^1000·10^46 >= 10^1000 with exact integers: true at exponent 46, false at 45.
- Claim 1: the same decoder is a straight-line program of 15 load/store pairs plus load, mask and store for the last word: 33 ops on every received word. 33^20 <= 2^101 < 33^20 · 2 (so L = 5 would be false). The experiment `ops_check.py` runs it on the reference counter for 22 received words (all report 33 ops and the correct output) and checks the inequality exactly.
- It is deliberately weak in error and cheap in decoding: one corner of the Pareto front. B displays as 1.0 (see section 3).

**Distinguishing cases (`samples/`):** each passes the structural preflight except the first, and each must be rejected. Each now carries both claims. In the first five, claim 1 is a plausible op count and the flaw is in claim 0. In the last four, claim 0 is the starter's valid error bound and the flaw is in claim 1.

| case | claim | intended outcome | what it catches |
|---|---|---|---|
| `below-converse` | B = 1/20 | refused by the score gate before review (`score_below_minimum`) | claims no code can achieve |
| `existence-only` | "there exists" a code with B = 585061/10^7, from the random-coding union bound | invalid | ensemble or existence argument instead of a fixed explicit code |
| `simulation-estimate` | an explicit (3,6) LDPC code from a seeded generator, BP decoding, B = 3/50 from 10^6 simulated frames | invalid | Monte Carlo confidence bound offered as the proof |
| `bit-error-rate` | repetition code, "error probability" 1/10 | invalid; the true block error is 1 − 0.9^1000 | bit-error rate passed off as block error |
| `checker-without-distance` | double-circulant code, bounded-distance decoding to radius 221, B = 281/5000; `bsc_exact` tail(222) passes; L = 1009 by exhaustive search | invalid | a protected-checker pass used without proving the minimum distance (sampled codeword weights are not a proof) |
| `expected-ops` | repetition code, first copy computed by a loop over disagreeing positions; L = 101/10 from the expected 180 disagreements | invalid; the worst case is 1000 disagreements | average-case cost offered as the worst case |
| `oversized-table` | repetition code, decoder is one lookup in a 2^1000-entry table; L = 6 | invalid | a table far beyond 2^32 words (and an address wider than a word) |
| `uncounted-loop` | repetition code, first copy by a 15-iteration loop; W = 33 counting only loads and stores | invalid; the loop's increments, tests and branches add 45 ops | loop control left out of the count |
| `decoder-mismatch` | repetition code; B for ML with random tie-breaking, ops counted for first-copy copying; L = 101/20 | invalid | an error bound and an op count for two different decoders |

The organizer list of what does not count is written as rules in `rules.md`, as the verifier's authoring guide asks. These cases are for builders to test against and are not part of the rules text.

**What could not be run here:**
- The reviewer's verdicts on the eight reviewed cases and the starter are intended outcomes, not observed ones: no live review has run.
- `test-checkers` and the sandboxed experiment step need Docker, which this build machine refuses. The checker's known answers were verified by running it directly, and its tail value by an independent floating-point sum (0.0561608773587…).

## Scope and provenance

This is a new challenge, not a reconstruction of a deployed one. The scientific contract is complete above; the files listed at the top implement it. The verifier pin is `da594aad5166918d4dd7d71630e1e649521e1953`, the challenge id is `bsc-2000-1000`, and the review standard is `bsc-2000-1000-review-v2`.

Intentional simplifications: decoding cost is counted in one fixed unit-cost model, so it ranks algorithms, not hardware. Encoding cost, memory use (beyond the 2^32-word cap), achieved error and code family are not ranked.
