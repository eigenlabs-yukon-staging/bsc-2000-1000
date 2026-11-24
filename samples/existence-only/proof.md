# A code with block error at most 0.0585061

Draw a codebook of 2^1000 codewords independently and uniformly from {0,1}^2000 and decode by maximum likelihood. The random-coding union (RCU) bound for this ensemble, computed exactly with integer and fraction arithmetic (as in the organizer's reference/bounds.py), is about 0.058506088975. The ensemble average of eps is therefore at most 0.0585061, so at least one codebook in the ensemble has eps <= 0.0585061. That code, with ML decoding, is the submission.


## Decoding complexity (claim 1)

ML decoding by exhaustive search over the 2^1000 codewords, stored as a preloaded table of 32 words per codeword. Per codeword: 32 loads of c, 32 loads of y, 32 xors, 32 popcounts, 31 adds, 1 compare, 2 selects, 1 counter add, 1 compare and 1 branch, at most 2^8 ops. The final 16 stores add 16. So W <= 2^1008 + 16 <= 2^1009.
