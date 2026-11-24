# Seven independently repeated message word planes

Write a 1000-bit message as a low-bit-first integer. The
first 1000 transmitted positions carry the message in order.
Set positions 1000 through 1023 to zero. For each word
j=0 through 6 of the message, copy its 64 bits to two further
transmitted words at indices 16+2*j and 17+2*j. Set all
remaining channel positions from 1920 through 1999 to zero.
The first copy makes the encoder injective. Exactly 448
message bits each have three disjoint transmitted positions.
The other 552 message bits each have one position.

The decoder copies the singly represented bits from the
first 1000 received positions. For each of the seven full
groups of 64 triply represented bits, it decides each bit
by majority. Three received words T,A,B provide the three
copies in the same bit order. Compute the 64 majority
decisions as

    C = T AND A
    D = T OR A
    D = D AND B
    T = C OR D

For each bit position, if T and A both set the bit then C
sets it. If exactly one of T,A sets the bit, D and then
the last OR set it exactly when B does. If T and A both
clear the bit then both intermediate terms clear it.
This proves the majority rule for every received word.

The counted program is the literal expansion of schedule()
in verify.py. The generator loops and decisions use literal
position indices to select which instructions to emit. The
emitted machine program uses one literal load for each of
the first sixteen received words, two more loads and four
logic operations for each of the first seven words, one
AND with 2^40-1 for the last output word, and sixteen
literal output stores. There are no runtime branches,
loops, address calculations or tables. Padding bits in the
last output word are always zero, irrespective of noise in
positions 1000 through 1023. The five registers T,A,B,C,D
are reused, all input addresses lie in 0 through 29, and
all output addresses lie in 32 through 47.
The worst-case count is exactly 16+7*6+1+16=75.
The exponent 3115/500, equivalently 623/100, satisfies
75^100<=2^623 by exact integer comparison.

For a triply represented bit, success means zero or one
of the three flips. Its probability is (9/10)^3 plus
3*(1/10)*(9/10)^2, which equals 243/250. For a singly
represented bit success means no flip and has probability
9/10. All these events concern disjoint channel uses, so
they are independent for each fixed message. The success
probability is constant in that message, hence the uniform
average error is exactly

    E = 1 - (243/250)^448 * (9/10)^552.

The declared calculation derives this exact fraction,
checks B=ceil(10^250*E)/10^250, checks
the exponent comparison, and exercises the schedule using
the reference word-RAM counter on sample inputs. The proof
gives the operation count and decoding rule universally.
