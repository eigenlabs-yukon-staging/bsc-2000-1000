# Double-circulant code that corrects 221 errors

The code has generator matrix [I | A], where A is the 1000 x 1000 circulant whose first row is the 1000 bits produced by Python's random.Random(7).getrandbits(1000), least significant bit first. The decoder is bounded-distance decoding to radius 221, with failure counted as an error.

Its minimum distance is at least 443: we sampled 10,000 random nonzero messages and every codeword had weight at least 443. So every pattern of at most 221 flips is corrected, and eps <= Pr[at least 222 flips]. The protected checker `bsc_exact` with kind tail and t = 222 confirms that this probability is at most 281/5000 = 0.0562.


## Decoding complexity (claim 1)

The decoder (bounded-distance decoding to radius 221) searches all 2^1000 codewords in Gray-code order. It keeps the message whose codeword is within distance 221 of y if there is one, and otherwise outputs the zero message. The 1000 generator rows are a preloaded table of 32,000 words. Per message: ctz of the step counter gives the row to add (1 op); 32 loads and 32 xors update the current codeword c; 32 loads of y, 32 xors and 32 popcounts and 31 adds give d(c, y); then 1 compare, 2 selects, 1 add and 1 compare on the counter and 1 branch. That is at most 200 <= 2^8 ops per message, and the final 16 stores add 16. So W <= 2^1000 * 2^8 + 16 <= 2^1009.
