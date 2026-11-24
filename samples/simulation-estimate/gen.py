"""Deterministic (3,6)-regular LDPC parity-check matrix by socket matching: n = 2000, 1000 checks, seed 1."""
import json
import random

n, m, wc, wr = 2000, 1000, 3, 6
sockets = [v for v in range(n) for _ in range(wc)]
random.Random(1).shuffle(sockets)
checks = [sorted(sockets[r * wr:(r + 1) * wr]) for r in range(m)]
print(json.dumps({"n": n, "checks": checks}))
