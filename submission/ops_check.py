"""Run the starter decoder on the bsc-ops-v1 reference counter (opmodel.py, a
byte-identical copy of the organizer's reference/opmodel.py) for several received
words, and report the op count, whether the output equals the first copy, and
the exact check W^2500 <= 2^12611 for the claimed bound 2^(12611/2500)."""
import hashlib
import json
import random
import sys
from pathlib import Path

# In the verifier sandbox the submission's files are at /snapshot/files/; locally, beside this file.
HERE = next(d for d in (Path("/snapshot/files"), Path(__file__).resolve().parent) if (d / "opmodel.py").is_file())
sys.path.insert(0, str(HERE))
from opmodel import Machine, INPUT_BASE, OUTPUT_BASE  # noqa: E402


def decode(m):
    # Straight-line program: no branches, so every received word costs the same.
    for j in range(15):                      # words 0..14 hold y_0..y_959 = first copy
        m.store(OUTPUT_BASE + j, m.load(INPUT_BASE + j))
    w = m.load(INPUT_BASE + 15)              # y_960..y_1023: y_960..y_999 are first copy
    m.store(OUTPUT_BASE + 15, m.and_(w, (1 << 40) - 1))


cfg = json.load(sys.stdin)
rng = random.Random(cfg["seed"])
words = [[0] * 2000, [1] * 2000] + [[rng.getrandbits(1) for _ in range(2000)] for _ in range(cfg["random_words"])]
counts, correct = [], True
for y in words:
    m = Machine(y)
    decode(m)
    counts.append(m.ops)
    correct &= m.output_bits() == y[:1000]
W = max(counts)
print(json.dumps({
    "opmodel_sha256": hashlib.sha256((HERE / "opmodel.py").read_bytes()).hexdigest(),
    "words_tested": len(words),
    "distinct_op_counts": sorted(set(counts)),
    "outputs_equal_first_copy": correct,
    "W": W,
    "W_pow_2500_le_2_pow_12611": W ** 2500 <= 2 ** 12611,
    "L_strictly_below_starter": 12611 * 20 < 101 * 2500,
    "W_greater_than_32": W > 32,
}))
