# Best [2000, 1000] code on BSC(0.1): error versus decoding cost

Carry 1000 bits in 2000 uses of a binary symmetric channel with crossover probability 1/10. Submit an explicit code (encoder and decoder) with two proofs: its average block error probability is at most B, and its decoder executes at most W <= 2^L operations on every received word, counted in the fixed bsc-ops-v1 model (64-bit word RAM, every primitive including loads, stores and branches costs 1). This is a Pareto challenge: both B and L are scored, and lower is better in each.

- No code can go below about 0.052548 in B (the sphere-packing converse); the verifier rejects any claim under 0.0525481183 before review. Random codes reach about 0.0585, but an existence argument names no code and does not count.
- L is between 4 and 4096. Exhaustive maximum-likelihood search costs more than 2^1000 operations; a first-copy decoder costs 33.
- The starting submission is a deliberately cheap, weak baseline: a repetition code with first-copy decoding, B = 1 - 10^-46 and L = 101/20 (33 operations).

Read rules.md. Edit only submission/: put two exact scored claims in submission.json (claim 0 gives B, claim 1 gives L, for the same code and decoder), the proofs in proof.md, and any code data or computation programs beside them. Declare every computation the proofs rely on as an experiment, which the verifier runs in a sandbox. `reference/opmodel.py` is the executable op counter; copy it byte-for-byte into submission/ to count your decoder on sample inputs. Run `llm-verifier check --root .` and `llm-verifier preflight --root . --submission submission` for structural validation. `python3 -I reference/bounds.py` reproduces the reference error values exactly.

The configured LLM reviewer judges the proofs; a person confirms official results. Local checks do not establish mathematical acceptance.

`samples/` holds test entries that must be rejected, one per likely misreading. BRIEF.md section 6 lists each with its intended outcome. `checkers/bsc_exact.py` is an optional protected exact-arithmetic tool, and `checker-tests.json` holds its known-answer cases.

<!-- llm-verifier:solver-guidance:begin -->
## Preparing a submission

Read the [solver guide](https://github.com/Layr-Labs/llm-verifier/blob/939a35dbb233c062ae514bd9abc62d3168eba1fc/docs/solving.md) and
[solver skill](https://github.com/Layr-Labs/llm-verifier/blob/939a35dbb233c062ae514bd9abc62d3168eba1fc/skills/llm-verifier-solver/SKILL.md)
from the same pinned verifier revision as this challenge.

- Use the challenge's verifier revision for `llm-verifier check --root .`.
- Run the final submission's declared computations through `intake`,
  `prepare-runtime`, and `experiments` without model credentials.
  `test-checkers` runs organizer fixtures, not your candidate's checks.
- Ask a fresh independent subagent to critique the complete final submission,
  protected rules, selected policies and evidence using the guide.
- Fix findings and rerun affected checks and review.
- Report results and omitted steps. Subagent critique is advisory;
  configured-provider review is separate. `simulate` does not judge correctness.
<!-- llm-verifier:solver-guidance:end -->
