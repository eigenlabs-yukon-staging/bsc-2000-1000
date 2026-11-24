# Four coordinate groups in 128- and 80-bit cyclic components

## Explicit encoder and information sets

Split a message into fifteen 64-bit components followed by forty bits.
Use length n=128 for the fifteen, and n=80 for the tail. There are
15*64+40=1000 information bits and 15*128+80=2000 channel uses.
All bits are numbered from zero, lowest first.

For p=n-1, the coordinates are F_p followed by infinity. Let G_p in
F2[X] be the monic gcd of X^p+1 with the sum of monomials whose
exponents are the distinct nonzero squares modulo p. Write k=n/2.
The k generator rows have bodies X^i G_p, i=0,...,k-1, and an
overall parity bit at infinity. The experiment computes their degree
and checks that elimination on the coordinate range k-1,...,p-1
gives its identity matrix. Denote this standard information set by I.
The encoder uses the I-systematic basis and sends its rows XORed
according to the message. Each block encoder is injective. Their
product therefore has exactly 2^1000 encoded words.

For the decoder, split the coordinates into four groups of size n/4.
Sort the coordinate indices by the SHA-256 byte strings of the
ASCII text four-qr:p:seed:j, with literal hyphens and colons as
in group_positions(), and unpadded decimal p, seed, j. Break any
hash ties by increasing j. Cut the result into four consecutive
groups and sort each group internally by coordinate index.
Seeds are 22 for p=127 and 14 for p=79. These fixed choices require
no probabilistic assumption. For each of the six pairs of groups,
form their union by concatenating the two groups in increasing
group order. The declared run eliminates the same generator rows
on each union and checks its identity restriction. This establishes
that every pair is an information set and fixes its systematic
basis.

## Exact spectra

The run checks that each generator row has weight divisible by four
and the rows are pairwise orthogonal. The rank has already been
checked by the identity on I, and k=n/2. Thus the code C equals its
dual. Every span word is doubly even by induction, since for two
orthogonal words the overlap cardinality is even and their XOR
weight is the sum of their weights minus twice the overlap.
The all-one word is orthogonal to all even words and belongs to C.

The run verifies primality and p=3 modulo four. It checks the
permutations x -> x+1, fixing infinity, and x -> -1/x, exchanging
zero and infinity. Apply each to each I-systematic row and
reconstruct the moved row from its I restriction. Equality means
both permutations preserve the code span. For each p, their
group contains every determinant-one projective permutation:
translation powers give U(a)=[[1,a],[0,1]], since the field is
prime, and the matrix S=[[0,-1],[1,0]] conjugates these to all
lower transvections L(a). Row additions and S reduce any determinant
one matrix to diag(a,1/a). These diagonals also lie in the group,
as U(a)L(-1/a)U(a)S=diag(a,1/a).

This group is transitive on unordered triples. An invertible
projective linear map sends any ordered triple to (infinity,0,1)
by sending the first two homogeneous lines to the axes and
rescaling the axes so that the third becomes 1. If its determinant
is nonsquare, compose with x -> 1/x, a permutation of the target
triple with nonsquare determinant -1. The latter is nonsquare
because the field multiplicative group has order p-1 not divisible
by four. The square determinant now allows rescaling the matrix
to determinant one. This proves the claimed triple action.
Transitivity on points and on pairs also follows directly from
the determinant-one matrices: inversion and translations move any
point to infinity, and translations fixing infinity move any
finite point to any other finite point.

Fix a support of weight w and let X be the size of its intersection,
after a uniform group action, with the standard information set I.
For j=1,2,3, transitivity gives
E[binom(X,j)] = binom(w,j)binom(k,j)/binom(n,j).
These factorial moments determine the first three ordinary moments.
They match those for a uniform w-subset of the n coordinates.
The latter has mean w/2 and third central moment zero by
interchanging the two halves of size k. Its variance is
w(n-w)/(4(n-1)), from summing the variances and pair covariances
of the indicators in sampling without replacement.

For weight four the variance is positive, so an image has X at
most one by mean two. For w=8,12 and (at n=128) w=16 the variance
is greater than one. Put Y=X-w/2. If all Y were at least -1, the
integer polynomial Y^3-Y would be nonnegative and vanish only
at -1,0,1. Its expectation is zero, forcing variance at most
one. Thus some image has X<=w/2-2. In particular any nonzero
word of weight below 20 for length 128 has an image with
1<=X<=6, and any nonzero word of weight below 16 for length 80
has an image with 1<=X<=4. The run enumerates the I-systematic
supports through degree seven and four respectively and checks
that every corresponding word meets the claimed distance.
The lower bound on X holds because I is an information set and
the group action preserves nonzero words.

For w=20 at length 128 the same method gives an orbit
representative with X<=7. If every X were at least 8, Y would
be at least -2 and the polynomial (Y+2)(Y-1)(Y-2) would be
nonnegative on those integers. Its expectation equals
4-Var(Y), strictly negative because 20*108/(4*127)>4.
This is a contradiction. Consequently the degree-seven search
finds a representative from every orbit of weight 20.
The run collects the weight-20 hits, applies the two verified
generators repeatedly with a queue and a set, and counts the
union. Positive powers of these finite-order permutations include
their inverses, so closure under the two gives precisely the
group orbits. All points visited are weight-20 codewords.
This computes A_20=56896 exactly for length 128.

The remaining spectrum is derived from exact duality equations.
For A_d denoting weight counts, define
K_w(d)=sum_j (-1)^j binom(d,j)binom(n-d,w-j).
The character sum over C at a word z equals 2^k if z is in
C^perp and zero otherwise. In the zero case, pairing with a
codeword having odd dot product cancels the sum. Sum over
weight-w words z and group the inner terms by overlap size with
a fixed weight-d codeword. By self-duality one obtains
2^k A_w = sum_d A_d K_w(d).
We know A_0=A_n=1, weights are divisible by four, and all
nonzero weights below the proved distance are absent. All-ones
complementation gives the mirrored coefficients. A_20 and
A_108 are also known for length 128 from the orbit computation.
The run uses exact rational elimination on the indicated
remaining columns. It verifies a pivot in each unknown and
consistent residual rows, giving the unique integral
nonnegative spectrum of each explicit code. No spectrum
coefficient comes from a numerical assumption or a citation.

## Decoder and success certificate

For a noise vector e, obtain its four group weights. For a
pair of groups with weights a,b, let u=min(a,b), v=max(a,b).
The pair is permitted if u is a valid zero-based index into
the profile and v <= profile[u]. The profiles for p=127 and
p=79 respectively are (4,4,4,3) and (5,4,3). A noise vector
is admitted if at least one pair is permitted. Both profiles
are nonincreasing in u, so if any pair is permitted then the
pair of two smallest group weights is permitted. Conversely
that pair is one of the six. This equivalence is used in the
volume calculation.

For each of the six information sets, compute the unique
codeword agreeing with the received word on that set. Add to
it every codeword corresponding to a permitted pair of mask
weights on its two groups. This produces exactly the
candidate codewords whose difference from the received word
is admitted via that pair: the noise on the information set
uniquely determines the difference from the agreeing codeword.
The zero mask ensures the candidate list is never empty.
Choose from all scans a codeword minimizing the weight of
its difference, with earlier scans and earlier positions within
a scan winning ties. Pair order and table order are fixed in
construction.py. Return its standard information set bits.

For any fixed received word y, its possible noise vectors
as the message varies are the coset y+C. The minimum admitted
weight depends only on this coset. With a uniform k-bit
message, the contribution to average success at y is
2^-k times the channel probability of the selected noise.
Independence of that probability from tie outcomes follows
because equal weights have equal channel probability.
Summing over the 2^k words of each coset gives that coset
one term 9^(n-t)/10^n, where t is its least admitted weight.

For each t let F_t be the number of cosets reaching by
weight t under the admitted noise restriction, and V_t the
volume of those admitted errors. The run computes V_t by
summing products of four binomial coefficients for the
four group weights, each over the n/4 positions of a group,
retaining only weights summing to at most t with permitted
two smallest values. For a nonzero codeword of weight d,
the number of ordered errors in the full radius-t ball whose
difference is that word is the intersection of the balls
about zero and that word. If j of its d positions are
selected and h outside are selected, both errors are in the
ball exactly when h<=t-max(j,d-j). This gives the
intersection(n,d,t) formula. Summing it with A_d for
0<d<=2t upper bounds the ordered distinct collisions among
the admitted errors. Call the upper bound M_t.

If m_s are coset occupancies among those errors,
sum_s m_s=V_t and sum_s m_s^2<=V_t+M_t. Cauchy-Schwarz
gives F_t>=ceil(V_t^2/(V_t+M_t)); the inequality
1>=m-binom(m,2) for positive integer m gives
F_t>=V_t-floor(M_t/2). For integer r>=1,
(m-r)(m-r-1)>=0 gives
r(r+1)*1_(m>0)>=(2r+1)m-m^2, including the zero case.
Summing yields F_t>=ceil((2rV_t-M_t)/(r(r+1))).
Use r=1+floor(M_t/V_t). The run maximizes these bounds and
the previous cumulative bound at each t, starting from
zero. Write the resulting cumulative lower bound as D_t.
It computes
A=sum_(t=0..tmax)(D_t-D_(t-1))*9^(n-t),
with tmax=22 and 18 respectively. This lower bounds the
success numerator because the coefficients of D_t in its
summation-by-parts expression are nonnegative, and replacing
D_t by F_t sums the success terms through tmax.
Later terms only increase success.

Uniform component messages and independent flips make the
sixteen component success events independent. If A and T
are the two computed numerators, full success is at least
A^15*T/10^2000. The stated rational error bound is checked
against that expression using exact integers.

## Word RAM program and count

The 15 main received blocks each occupy two entire words,
at literal addresses 2i,2i+1 for i=0,...,14. The last one
occupies words 30 and 31, with word 31 zero above its low
16 bits. At each block, load the two words Y0,Y1.
The output blocks are one word each, with the tail occupying
the low 40 bits of word 47. No crossing is needed.

Starting at address 48, concatenate the six delta tables
for the length-128 code and then the six for the tail.
For each information set, its table has two words per
permitted mask. Store the corresponding systematic-basis
codeword low word first, followed by its remaining high
bits. The masks are enumerated in the order of increasing
first weight, then increasing second weight, then increasing
lexicographic support tuples for the first group and second
group. The exact count per table is
N=sum_(permitted a,b) binom(n/4,a)binom(n/4,b), as in
table_entries() with delta_words() giving every entry.
Each block of the same length uses the same six tables.
They are fixed offline. The sum of the lengths and the
48 original words is 827218848, within the model.

This paragraph and the SCAN listing specify the whole
program. All block and information-set loops and loops over
basis coordinates are code-generation loops, emitting repeated
literal-address operations. The only runtime loops are the SCAN
loops below. Moves and assigning constants are free as specified
in the protected model. At the start of a block, move n+1
into BestW, and zero into Best0,Best1. For each of the six
information sets, move zero into B0,B1. For each fixed
coordinate in its ordered pair of groups, shift its received
word Y0 or Y1 right by the coordinate's location in that word
and AND with 1. Use that bit with two selects to choose the
fixed corresponding basis row's low and high words or zero.
XOR these into B0,B1. These are six primitives per coordinate
and exactly form the agreeing codeword. Compute E0 and E1 by
XORing the received words with B0 and B1. Move the table's
literal base into Addr and execute this loop:

    SCAN:
        D0 := load(Addr)
        Next := add(Addr,1)
        D1 := load(Next)
        N0 := xor(E0,D0)
        N1 := xor(E1,D1)
        C0 := popcount(N0)
        C1 := popcount(N1)
        Wt := add(C0,C1)
        Take := ltu(Wt,BestW)
        BestW := select(Take,Wt,BestW)
        Best0 := select(Take,N0,Best0)
        Best1 := select(Take,N1,Best1)
        Addr := add(Addr,2)
        More := ltu(Addr,End)
        branch if More to SCAN

The not-taken branch falls through to the next emitted
instructions. End is the table's literal first address
after its last pair. The loop executes exactly N iterations
at 15 primitives each and no other test or jump. The two
noise words and weight update exactly when the new candidate
strictly improves the previous best. Across the six scans
this implements the tie-broken decision above.

XOR Y0,Y1 with Best0,Best1. Shift the resulting low word
right by k-1 and its high word left by 65-k, then OR them.
This extracts the I bits into the low k bits of D. For k=40
AND with 2^40-1 to enforce padding. For k=64 the low word
already holds exactly the desired bits. Store D at literal
address 32+i for component i. Halt after the last store.
There are at most 25 data registers including all named
temporaries, with no working memory beyond the given ranges.

Thus the exact count for an n-bit block is
2+6*(6k+2+15N)+5+1_(k<64)+1. All 16 blocks together
execute W=92104426581 primitives for every input. The
declared integer power comparison is W^1000 <= 2^36423,
so L=36423/1000. The all-input count follows from the
literal schedule and its bounded scans. The declared
experiment verifies the numbers and does not rely on
attempting a multi-billion-instruction sample decode.
