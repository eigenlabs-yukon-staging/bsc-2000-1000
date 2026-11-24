# Information charts with interpolation multiplication

## Encoder, decoder and chart family

Let F be the field of size 2^50 with its polynomial
basis over P=X^50+X^9+X^2+X+1. The certificate verifies
X^(2^50)=X modulo P and gcd one for X^(2^25)+X and
X^(2^10)+X with P. To see why these checks suffice,
an irreducible factor of degree j has a quotient field
where Frobenius has order j. The jth power fixes all
elements by the multiplicative group order. A smaller
positive power cannot fix all elements by the root
bound. Fixing X fixes that entire quotient because
X generates it. Therefore the first check forces
every factor degree to divide 50, and the two gcds
exclude all proper such degrees. This proves P
irreducible. The same argument over F proves
irreducibility of

    R=T^10+T^3+108*T+184

by checking the q0^10 power of T modulo R, where
q0=2^50, and the monic gcd with T^(q0^j)+T for
j=2 and 5. The constants 108 and 184 mean their
binary polynomial basis representations in F.
Polynomial Euclidean division in verify.py uses the
proved arithmetic on F. Thus K=F[T]/R is a field of
size q=2^500. Write theta for the residue of T.

Represent each element as its ten 50-bit
coefficients on 1,theta,...,theta^9, low coefficient
first. Each coefficient is written in low-to-high
order in the F polynomial basis. Fix lambda=7 in F,
and h=(q-1)/3. The certificate checks that
zeta=lambda^h is a nontrivial cube root of unity,
so lambda cannot be a cube.

From two field elements a,b in the 1000-bit input
the encoder sends four consecutive full groups

    (a,b,c,d),   c=a*b,   d=a^3+lambda*b^3.

It is systematic on all 2^1000 messages. The
decoder considers only messages with a and b
nonzero. It chooses an eligible candidate of minimum
full Hamming distance, with ties to the smallest
integer a+2^500*b in the given binary order.
Eligibility is the union of the following charts.

- Six ordinary pairs (a,b),(a,c),(a,d),(b,c),(b,d),(c,d)
  take their two whole groups, with total distance
  caps (99,99,101,99,101,101).
- Each split chart takes either full systematic group
  v=a or v=b, any five of the ten coefficients of
  the other systematic group u, and any five coefficients
  of c. The total known bits are 1000. The chart cap
  is 99. It also requires that the recovery submatrix
  for that v and these index sets be invertible.

There are 2*252^2=127008 split charts and six
ordinary charts. The recovery condition is defined
below and depends only on the outer field and the
coefficient indices, all of which each key supplies.
The eligible set is nonempty since the two received
systematic groups can be made nonzero with at most
two bit changes and then qualify via (a,b).
The analysis uses a subset of messages with every
recovery submatrix invertible.

## Parity differentials and multiplicities

For a fixed nonzero systematic difference (s,t), the first
parity difference is a*t+b*s+s*t. The second is

    s*a^2+s^2*a+s^3+lambda*(t*b^2+t^2*b+t^3).

These are affine F2-linear functions of a,b. Suppose s!=0.
Fixing the first parity difference leaves b affine in a with
its linear coefficient t/s. Substituting into the second
gives an affine function of a whose linear part is
A*(a^2+s*a), where A=(s^3+lambda*t^3)/s^2.
This scalar is nonzero, since if t!=0 an equality to zero
would make lambda a cube, and if t=0 it equals s.
The kernel is {0,s}. If s=0 then t!=0. Fixing the first
parity difference determines a, while b maps into the
second via lambda*(t*b^2+t^2*b) with kernel {0,t}.
So in every direction each occupied parity difference
has exactly two preimages. The parity difference for
a uniform message is uniform on an affine hyperplane
in the 1000-dimensional binary parity space.

We bound how often different directions can have the same
binary mask for this hyperplane. The F2 trace on K is
Tr(x)=sum_(i=0..499) x^(2^i). It is F2-linear with values
in F2, and nonzero by the root bound. Its coordinate pairing
is nonsingular. If Tr(z*x) were zero for all x for z!=0,
the trace would vanish everywhere since multiplication by z
is surjective. Thus every binary parity mask has a unique
trace pair (alpha,beta) such that its form on (c,d) is
Tr(alpha*c+beta*d).

The mask annihilates the linear part of the derivative in
direction (s,t) precisely when the two coefficients of
a and b, after absorbing squares under the trace, vanish:

    alpha*t + beta*s^2 + sqrt(beta*s) = 0,
    alpha*s + beta*lambda*t^2 + sqrt(beta*lambda*t) = 0.

Here the field square root is unique, and the identity
Tr(z*a^2)=Tr(sqrt(z)*a) follows by squaring under the trace.
If beta=0 and the mask is nonzero, these equations force the
zero direction. If alpha=0 and beta!=0, squaring the
equations gives one nonzero quartic each in s and t,
so there are at most sixteen pairs, including zero.
When alpha and beta are nonzero, square the equations.
The first expresses t^2 as A0*s^4+B0*s with A0!=0.
Substituting its square into the second expresses t as
A1*s^8+B1*s^2 with A1!=0. Squaring and comparing with
the first expression for t^2 gives a nonzero degree
sixteen polynomial in s, and t is uniquely fixed by s.
Thus the total solution count is at most sixteen including
zero. A given nonzero mask can therefore occur for at
most fifteen nonzero directions. The hyperplane for a
given direction has a unique nonzero mask because of its
proved rank.

## Good messages and exact admitted volume

For field v, its multiplication matrix M_v over F
has column j given by the coefficients of theta^j*v.
Every entry is a linear form in the ten coefficients
of v. Consider all 63504 choices of five rows and
five columns. The determinant of each selected
square submatrix is a polynomial of degree at most
five in those ten coefficients.

The certificate shows all these determinant
polynomials are nonzero by evaluating them at the
single concrete WITNESS given in construction.py.
It computes each k-by-k minor for k=1..5 from
the preceding size by expansion along the largest
row index. In characteristic two every sign in
that expansion is plus. Multiplication by each
of the witness matrix's fixed entries is factored
into seven byte tables by linearity. The program
computes these tables from the base-field
polynomial and counts nonzero minors at size five.
The exact count is 63504.

A nonzero polynomial of degree at most d in n
variables over a field of size U has at most
d*U^(n-1) zeros. The univariate root bound starts
the induction. View the polynomial as one in the
last variable of degree k whose leading coefficient
has degree at most d-k. When that coefficient is
nonzero, there are at most k choices of the last
variable that are roots. On its zero set, induction
bounds the number of prefixes by (d-k)*U^(n-2)
and at most U final values per prefix. This proves
the bound. Call v analysis-good if it has all the
stated minors nonzero. Since each zero set has at
most 5*2^450 elements, there are at least
q-317520*2^450 analysis-good values, all nonzero.
For Q=2^1000 set G=(q-317520*2^450)^2. The true
number Gtrue of messages with both fields
analysis-good is at least G.

For such a message, chart admission depends only
on the noise vector. Write E_t for the admitted
errors of total weight at most t<=240, and V_t
for their number. For a received y let m_y be
the number of analysis-good eligible messages
within distance t. Then sum_y m_y=Gtrue*V_t.

Let A,B,C,D denote the noise weights in the four
whole groups. For X any group let v_X be the sum
of the five smallest chunk weights of that group.
Write s=A+B, p=C+D, m=min(A,B) and
u=min(A+v_B,B+v_A). Then all charts fail exactly
when

    s>99, p>101, m+C>99, m+D>101, u+v_C>99.

Indeed the optimal five fixed-size chunks from
each of the two chosen groups have exactly that
minimum sum, and any two choices of index sets
are included. For total weight through 240,
the failures satisfy 100<=s<=138 and
102<=p<=240-s.

Define T(w,k) as the number of masks on a single
whole group with total weight w and sum of five
lightest chunks at least k. It is binom(500,w)
at k<=0 and zero for k>floor(w/2). To compute
it, let r be the fifth smallest chunk weight.
If exactly j chunks are strictly below r and e
equal r, then j<=4, j+e>=5 and 10-j-e chunks
are strictly above r. For fixed r,j,e, choose
these chunk positions by the multiplier
binom(10,j)*binom(10-j,e). For the j low chunks
the polynomial

    (sum_(i=0..r-1) binom(50,i)*z^i)^j

gives counts by their total weight wlow. Likewise
for the high chunks use the polynomial containing
i>r to degree 50 raised to the number of high
chunks. Their total weight is whigh. The equal
chunks contribute binom(50,r)^e choices. For
these terms the whole weight is w=wlow+e*r+whigh
and the five-smallest sum is wlow+(5-j)*r.
There are at least six chunks of weight at least r,
so for any whole weight w<=140 only r<=floor(140/6)
matters. All relevant whole weights for counting
rejection have that bound by the s,p ranges above.
The declared calculation accumulates this histogram
and verifies that its row mass equals binom(500,w).
A suffix sum takes the count T(w,k).

For fixed A<=B, the count of masks with u>=j is
T(A,j-B)*T(B,j-A). Differences between consecutive
thresholds give the equality counts. The range
A..min(A+floor(B/2),B+floor(A/2)) contains all u.
Double the count when A<B to include swapped
groups. Holding A,B,u and total p fixed, the
number of remaining failures is

    sum_C T(C,100-u)*binom(500,p-C),

over C>=0, p-C>=0, C>99-m and p-C>101-m.
Combining these counts with the ranges of s,p
gives exact rejected shells. Subtract them from
binom(2000,t) for each total weight t and cumulate
to find V_t.

To bound collisions, for each error in E_t allow the
competitor to have any noise of total weight at most t.
This only enlarges the number of ordered distinct pairs.
For a fixed first error and fixed systematic competitor
noise of weight j, the parity noise lies in an ordinary
1000-bit ball of radius ell=t-j. Let J be the volume of
this ball. Write D0=42. A binary parity mask is balanced
for this bound when at least one of its two 500-bit parts
has at least D0 ones and D0 zeros. Such a mask has at
least D0 ones and zeros across the parity ball. The
signed shell count of its mask of weight w is the
coefficient of z^k in (1-z)^w(1+z)^(1000-w).
Pair D0 factors of each sign. The resulting signed
coefficients are bounded in absolute value by the
coefficients of (1+z^2)^D0*(1+z)^(1000-2*D0) by
the triangle inequality. Thus the absolute signed
ball count is at most H, the sum of those coefficient
bounds through ell. Translation by the fixed first
parity error multiplies the signed count by a sign.
The affine hyperplane selected by the direction
meets the translated ball in at most (J+H)/2 values.
For an unbalanced mask use the trivial bound J.
The number of such directions is at most

    Kbad=15*(2*sum_(i=0..D0-1) binom(500,i))^2

by the multiplicity proof in either chosen binary basis.

The systematic noise of a competitor determines the
direction when the first error is fixed, by binary XOR.
For each weight class j there are N=binom(1000,j) such
patterns, with at most min(N,Kbad) unbalanced directions.
The zero direction is excluded but can increase the
upper count if included in N. From the mass 2/Q on the
affine hyperplane, the expected number of competitors
over all uniform first messages is at most A_t/Q,
where A_t sums N*(J+H)+min(N,Kbad)*(J-H) over j=0..t.
The signed coefficient bounds satisfy H<=J because
their polynomial replaces (1+z)^(2*D0) factors by
(1+z^2)^D0, which is coefficientwise no larger.
Restricting either message to the good set can only
decrease the count. Thus sum_y m_y^2<=V_t*(Gtrue+A_t),
where Gtrue is the actual number of analysis-good messages.

Let F_t=|{y:m_y>0}|/Q. By Cauchy-Schwarz, use Gtrue first
to give F_t>=V_t*Gtrue^2/(Q*(Gtrue+A_t)). For x>=0,
x^2/(x+A_t) is nondecreasing, so Gtrue can be replaced
by the lower bound G. For integer l>=1, the inequality
(x-l)(x-l-1)>=0 for positive integer x gives
l*(l+1)>=(2*l+1)*x-x^2. It also holds with zero
left side at x=0. Summing this form over y and using
both moment estimates gives
F_t>=V_t*(2*l*G-A_t)/(Q*l*(l+1)), since replacing Gtrue
by G weakens its numerator.
The exact calculation uses l=1+floor(A_t/G), floors
both bounds, and takes their maximum with the previous
lower bound. Let D_t be this nondecreasing sequence,
with D_(-1)=0.

The cumulative number of received words with winning
distance at most t, normalized by Q, is at least F_t.
This follows because the decoder considers every
analysis-good eligible codeword and chooses a minimum. Each
received word with winning distance t contributes
Q^-1*9^(2000-t)/10^2000 to average success, because its
selected message has exactly one corresponding noise
pattern including on tied distances. For t=0..240 the
success lower bound is the sum of (D_t-D_(t-1))*9^(2000-t)
divided by 10^2000. Summation by parts gives nonnegative
coefficients for each cumulative D_t, so the actual
cumulative counts can only increase this expression.
The certificate finds one minus this lower bound at most
55632661054900983/1000000000000000000. This normalization
includes the entire uniform message distribution.

## Solving the six ordinary charts

The chart for a specified pair of group values produces
exactly its compatible (a,b) with a,b nonzero. Chart pairs
are visited in lexicographic pair order. On each chart
enumerate all two-group received perturbations of weight
at most that pair's cap. Retain only candidates with
both systematic fields nonzero. Duplicate messages do not alter
the minimum or numerical tie rule.

Pairs (a,b), (a,c) and (b,c) determine the message directly:
skip keys with a required nonzero value zero, and for
(a,c) or (b,c) recover the other component by division.
For (a,d) find all roots of b^3=lambda^-1*(d+a^3).
For (b,d) find all roots of a^3=d+lambda*b^3.
Skip a zero known systematic group or zero target in these
two cases.

For (c,d), skip c=0 or d=0. A good message has both
c!=0 and d!=0, because a^3+lambda*b^3=0 would make lambda
a cube. Define k=lambda*c^3. The two cubic quantities
a^3 and lambda*b^3 must be the roots of
A^2+d*A+k=0, since their sum is d and their product is k.
Compute a root w of w^2+w=k/d^2, or return no candidate
if none exists. Put A=d*w. It is nonzero and the other
root is A+d. One need only raise this one A to the cube
root exponent e below. Write t=A^e. If t^3=A, use
a=t, b=c/a. If t^3=zeta^2*A, then
(A^h)^2=zeta^2, so A^h=zeta. In this case let
b=t*lambda^(-e), whose cube is
A*zeta^2/(lambda*zeta^2)=A/lambda, and set a=c/b.
In the remaining case there is no good solution. Indeed,
for any good solution A must be either the cube a^3 or
lambda times the cube b^3, and those two cases force the
two stated checks on (A^e)^3. Conversely, in the first
case lambda*b^3=k/A=A+d, and in the second
a^3=k/A=A+d, so the equations hold. The cases are
disjoint since lambda is a noncube. Rotate the found base
(a,b) by (zeta*a,zeta^2*b) and (zeta^2*a,zeta*b).
These three pairs exhaust all roots of the chosen cubic
and their corresponding division values.

The macro for cube roots uses e=(2h+1)/3, computed exactly.
For a cubic z!=0, z^h=1 and (z^e)^3=z. The same identity
holds at zero. A proposed root x=z^e is checked by cubing.
If valid, its multiples by zeta and zeta^2 give all
roots by degree, with duplicates only at zero. In the
charts admitting cube roots the target is nonzero.
For its fixed addition chain use
e=1+56*sum_(i=0..82) 2^(6i).
Compute z^7 as z^2*z^4*z, raise to 2^3 to get z^56,
then take its geometric power of length 83 with step 6
and multiply by z. For inversion of x!=0 use
x^-1=(x^(2^499-1))^2. This takes geometric length 499
with step one, then one squaring. To compute a geometric
power of length L of input T, set G=T and a=1.
For each remaining binary digit of L after its leading
one, replace G by G*G^(2^(step*a)). If the digit is one,
also multiply by T^(2^(2*step*a)). Update the compile-time
a to 2a plus the digit. The identities follow by splitting
the geometric exponent into adjacent intervals.

The linear map taking z to a possible solution of
w^2+w=z is AS_map(). The code computes a basis of the
image of x^2+x with tracked preimages, using descending
pivot positions. Reducing z against the fixed pivots
and accumulating their corresponding preimages is linear
over F2. On the image the residual is zero, and its
output is a preimage. The program verifies the equation
after this map, so for z outside the image it discards
its result. The two preimages for an attained z differ
by one, since the kernel is {0,1}. The declared computation
checks the AS map on the image basis, the exponent
identities and the chart rules on sample messages.

## Solving a split chart

Let v be the known full field and u the other.
The inner key supplies five coefficients P of u
and five coefficients Z of c=v*u, each in
increasing index order. Let U be the other
five coefficients of u in index order. With
rows and columns selected from M_v as specified,
the key implies

    M*U + J*P = Z.

Rows of both matrices select the c coefficients.
Columns of M select unknown u coefficients,
while J uses the known ones. If M is singular,
skip this outer key for the chart. Otherwise
there is exactly one U=M^-1*J*P+M^-1*Z.
The full v and u determine c and d. The last
group equals Dbase+alpha*u^3 with Dbase=v^3
and alpha=lambda if v=a, and Dbase=lambda*v^3
and alpha=1 if v=b. This solves exactly the
split chart for a,b nonzero. Any eligible
message supplies a possible key in its chart,
and every retained codeword matches all known
key bits. Together with the ordinary pairs,
minimum-distance scoring then gives the
analyzed decoder.

## RAM layout and field primitives

Each K value uses ten words, one per 50-bit
coefficient, with the top fourteen bits of every
word zero. This is padded working storage. Input
and output pack the chunks back-to-back.
ARRAY_LIST defines 67 arrays of ten words from 48.
The four traversal stacks NEXT,LAST,INEXT,ILAST
each use 102 words at 718,820,922,1024. MULT has
sixteen single-word rows at 1126. The two fixed
reduction tables have two and sixteen entries at
1142 and 1144. MAP_BASE is 1160 and holds the 41
maps of map_tags, with seventy byte slices of
256 entries of ten words each.

The fixed chart descriptors start at CHART_META
at 7348360. Each row has 16 words as in descriptor,
the orientation, known systematic index set,
missing index set and selected product index set.
The 127008 rows are in SPLIT_CHARTS order.
DIRECTIONS starts at 9380488 with 500 rows of
twenty words. SAVED at 9390488 has 102 rows of
thirty words. AFFINE at 9393548 has 500 slots
of 2662400 words, each with 40 variable-size
slices specified below. These are writable scratch.

EVAL_A, EVAL_B and EVAL_PRODUCT are three 19-word
scratch arrays at 1340593548, 1340593567 and
1340593586. Fixed interpolation-multiplication
tables EVAL_TABLE and INTERP_TABLE start at
1340593605 and 1739247045, respectively. Their
layout and contents are given below. The exclusive
extent is 2137900485, within 2^32, and all ranges
are disjoint.

The macro loops over literal words, fixed choices
and power expansions below run at code-generation
time unless a runtime counter is specified.
Array references in macros denote word addresses.
A K COPY costs 20 (one load and one store per
word), and a K XOR costs 40 (two loads, XOR and
store per word). A field ZERO costs 21, using ten
loads, ten ORs into a zero accumulator and an
equality. EQ costs 41, using two loads, XOR and
OR per word and one equality. A following branch
is charged separately.

For a one-word F element w use BSHIFT_s for
fixed s=1 or 4. Emit

    c := shr(w,50-s)
    w := shl(w,s)
    w := and(w,2^50-1)
    addr := add(REDSHIFT_s,c)
    h := load(addr)
    w := xor(w,h)

Here the table entry h is the carryless product
of c and P+X^50, fitting in a word and implementing
the modular reduction. Thus BSHIFT costs six
operations and multiplies by X^s in F.
BMUL on F loads its first operand, stores zero
and the operand in rows 0 and 1 of MULT, and
applies BSHIFT_1 three times with each result
stored in rows 2,4,8. For the remaining eleven
nonzero rows in increasing order, compute the
XOR of the rows for the least set bit and the
remaining bits, with two loads, XOR and store.
This defines all nibble multiples of the first operand.
Start an accumulator at zero and for each of the
13 literal nibbles of the second operand, high
to low, load, shift, mask, multiply by one and
add MULT to obtain a row address. Apply BSHIFT_4
to the accumulator and use add of literal zero
to the row address, load and XOR. These operations
are Horner's rule. Store the resulting F word after
all nibbles. The count is
1+2+3*(6+1)+11*4+13*(5+6+3)+1=251.
Because the final store follows all operand reads,
the result may alias an operand. MULT is disjoint.

A fixed linear APPLY on K uses ten accumulators
initially zero. For each of the seventy literal
bytes, read its coefficient word, shift and mask,
multiply by ten and add the fixed slice base.
For each of ten word indices add it to that
address, load and XOR into its accumulator.
Store the ten accumulated words at the end.
The last byte of each coefficient has two live
bits, with the rest padded by zero. Each entry
is the map applied to the selected polynomial
basis coefficients for that byte. The functions
map_columns and map_words give the entire table.
The count is 70*(5+30)+10=2460. The
instructions above are all scalar.
Input and output may alias after the final read.

To multiply two K elements, first evaluate each
coefficient polynomial at the 19 distinct points
rho=0 through 18 in F, with these literal binary
representations. The degree of their product is
at most 18, so 19 point values determine it.
For the point products use nineteen BMUL, one
for each point. Interpolate by the Lagrange basis
polynomials

    L_i(T) = product_(j!=i) (T+rho_j)/(rho_i+rho_j).

They satisfy L_i(rho_j)=1 when i=j and zero otherwise.
Because the difference of two polynomials of degree
at most 18 with nineteen roots is zero, summing point
products times L_i gives exactly the full product.
Reduce the result modulo R. The reduction is linear
over F, so the interpolation map can use the fixed
vectors of coefficients of L_i modulo R. The function
interpolation_columns computes these vectors.
The fixed evaluation map uses the respective powers
rho_i^j in evaluation_columns. Both maps are linear
over F as well as over F2.

For each input coefficient, use three bit slices
of widths 20,20,10 at positions 0,20,40.
For one coefficient at a time reserve consecutive
slices in that order, with 2^bits entries per slice.
Each EVAL_TABLE entry has 19 words and gives the
evaluation map's output contribution from that
coefficient slice, which is the coefficient value
shifted into its actual bit position times each
rho_i^j. Each INTERP_TABLE entry has ten words
and gives the interpolation map's contribution from
that single point-value slice, using its corresponding
column. linear_slice_entry specifies both. These
ranges each contain 19*10*(2^20+2^20+2^10) words.

To apply either map, initialize its accumulators
to zero and loop over literal slices. Per slice
load its one input coefficient, shift and mask,
multiply by the output row length and add its
literal slice base. For each output index add
it to the address, load and XOR into the
accumulator. Store the accumulated words.
Evaluation therefore costs 30*(5+57)+19=1879
for either operand. Interpolation costs
57*(5+30)+10=2005. Emit these two evaluations
into EVAL_A and EVAL_B, nineteen BMUL products
into EVAL_PRODUCT, then the interpolation.
The full count is 2*1879+19*251+2005=10532.
We use the uniform MUL cap 11000. Its scratch
is disjoint from operands, which are read before
any result writes, allowing aliased outputs.
The maximum live map uses 19 accumulators.

CUBE uses APPLY of squaring into X2 followed by
MUL with its input, at most 13460. Use X4 for
scratch when required for alias separation.
INV uses one COPY, 14 MUL and 15 APPLY, costing
at most 190920. ROOT_POWER uses one COPY, 13 MUL
and 13 APPLY including its cube in CHECK, at most
175000. ROOT adds EQ for a total at most 175041.
These implement the geometric and root exponent
chains above. Use X2,X4,T for the initial cubic-root
powers and G,H,J for their geometric chains, preserving
other chart arrays throughout. Every count includes
address arithmetic and all field data movement.

At each candidate score, test both systematic
values for nonzero using ZERO and branches at
cost 44. Across the forty field words load each
candidate and corresponding received word, XOR,
popcount and add to a running distance, for 200
operations. Compare the twenty systematic words
in low-to-high order against BEST, emitting one
load of the stored word, ltu, eq and select to
replace the strict-smaller flag only for a differing
word. This costs 80 and realizes the tie rule
because padded fields preserve numerical order.
Compare the full distance and combine the
tie flag with two comparisons, AND, OR and select
for the new best distance. At each best word load
both prior best and candidate, select and store
for another 80 primitives. Total SCORE is 365.
The combined guard and score bound used below
is 425.

## Ordinary chart bodies

Each ordinary chart has an outer key of one full
group, the first named except for (c,d), whose
outer key is d. Copy KEY0 to its V array and branch
on zero to outer extension. Once per outer visit
compute KNOWN=a^3 for (a,b) and (a,d). For (a,c)
also compute a^-1 into INV. For (b,c) use KNOWN
=lambda*b^3 and INV=b^-1. For (b,d) compute only
KNOWN=lambda*b^3. In (c,d) compute d^-1 into INV
and its square into INV2. Including the guards
and copies these consume fewer than 600000
operations per outer visit.

At each inner visit COPY KEY1 into its V array
at cost 20. For (a,b), guard b, compute c by MUL
and d by CUBE(b), APPLY of lambda and XOR with
KNOWN. For (a,c) and (b,c), guard c, recover the
missing systematic group by MUL(INV,c) and form
d by CUBE of it, with lambda map only for b, and
XOR with KNOWN. Score. These first three body
types cost at most
20+22+11000+13460+2460+40+425=27427,
so the cap is 32000.

For (a,d), XOR d and KNOWN into TARGET and apply
the lambda^-1 map in place. For (b,d) omit that
map. Guard zero TARGET, then use ROOT with
failure branch. COPY the base root and use
APPLY of zeta and its square for the two rotations
into the missing systematic group. In each case
compute c by MUL and score. The bound is
20+40+2460+22+175041+1+20+2*2460+3*(11000+425)
=216799, below 220000.

For (c,d), guard c nonzero, and form K=lambda*c^3
using intermediate C1. Multiply K by INV2 into
TARGET. Apply AS to TARGET into ASROOT, square
ASROOT into TEMP and XOR ASROOT into TEMP. Test
EQ with TARGET and branch to extension on failure.
MUL(d,ASROOT,A0) forms A. Use ROOT_POWER(A0,ROOT).
Use EQ(CHECK,A0) and branch to ordinary orientation
if true. Otherwise use APPLY of zeta^2 to A0 into
TEMP, EQ(CHECK,TEMP) and branch on failure.
On success compute BASEB by APPLY of lambda^-e
to ROOT, then INV(BASEB,INVA) and MUL(c,INVA,BASEA).
Jump to shared scoring. In the ordinary orientation,
COPY ROOT into BASEA, then INV(BASEA,INVA) and
MUL(c,INVA,BASEB). Score by two COPY into V0,V1
for the base pair, then four APPLY of zeta or
zeta^2 for the pairs of rotated fields. Each pair
is scored separately and c,d remain unchanged.
The larger branch path, using the extra
orientation tests, map and jump, costs at most

    20+22+15920+11000+2460+2460+40+42+11000
    +175000+42+2460+42+2460+190920+11000+1
    +40+4*2460+3*425 = 436044.

The cap is 440000. Candidate guards branch past
the corresponding evaluation on failure, leaving
the remaining enumeration and scores intact.

## Split chart setup and table generation

A chart descriptor supplies the five known
systematic coefficient indices, their five
complementary indices and the five selected c
indices. The outer key KEY0 is v. At its outer
visit copy to the appropriate V field and guard
zero. The column arrays MATRIX0..MATRIX9 are
theta^j*v, computed by COPY for j=0 and fixed
APPLY maps otherwise. Each selected coefficient
of M or J is accessed by adding the row index to
the ten-word column address. The latter is
formed by multiply by ten and add of MATRIX0
from the selected column index. The indices
are read from their fixed descriptor offsets.

Use GAUSS0..GAUSS4 to hold the five rows of
the 5-by-10 augmented matrix [M|I]. For each
literal pivot column j=0..4, inspect the entries
j through 4 in that column in order, using load,
equality with zero and a branch for each.
If all vanish, skip the outer key. Otherwise
swap the found row into row j by ten pairs
of loads and stores. Distinct rows are in the
ten-word array region and the pivot row number
can be used with multiply and add for its base.
Copy the pivot entry to a fixed temporary and
use the fixed binary chain for exponent 2^50-2,
at most 100 BMUL operations, to compute its
inverse in BASEINV. Multiply each entry in
row j by this inverse by BMUL. For each
other literal row, copy its column-j entry
to a scratch word and replace its ten entries
by their XOR with that entry times the
corresponding normalized row-j entry. The
usual row-operation identity leaves the
augmented matrix [I|M^-1] exactly when
M is invertible. Copy columns of M^-1
into the first five words of COEF5..COEF9,
and compute columns of M^-1*J into the first
five words of COEF0..COEF4 by sums of five
BMUL products per entry. Each COEF row is
thus the missing-coefficient change for its
particular key lane. Use the first word of
ROW_TMP or BASEINV for temporary base products
as permitted by data lifetimes.

For the inner root read the five known
received coefficients of u and the five
selected received coefficients of c, using
the descriptor and the fixed Y arrays. Store
them in order in KEY1, and copy the first five
to the corresponding V field positions. The
remaining five positions are dot products of
the ten KEY1 words with their COEF entries.
Compute them by BMUL and XOR. Compute c=MUL(v,u)
and d by two cubes, lambda map on the appropriate
term and XOR. This initializes all four V
arrays for the actual inner root. The inverse,
coefficient calculation and root setup use
fewer than 1000 BMUL, 20 APPLY, two CUBE,
one MUL and a fixed number below 5000 of
other scalar operations. For example the
pivot phases use at most 500 inversion BMUL,
50 row scaling and 200 elimination; the
coefficient products add 125 and the root
products add 50. All guards, address operations
and copies in that setup and the initial matrix
columns fit below 1500000 operations.

Generate a direction row for each possible
inner bit i. It changes the key chunk
j=floor(i/50) by f=X^(i mod 50). A set-bit
shift from literal 1 creates f. For the missing
five coefficients of u, use BMUL of f with
the five entries of COEF[j]. Set the five known
coefficients to zero except for the appropriate
f when j<5. Place those literal or addressed
words into DELTA to obtain e, the full change
of u. They form a linear change because
U=M^-1*J*P+M^-1*Z. Use MUL(v,e) for the
change of c. The row at DIRECTIONS+20*i holds
the two ten-word arrays in that order. The
runtime emits ten literal lane loops each of
50 values, choosing the COEF array and which
known index to update by literal lane choice.
Its initializations, comparisons, branches,
counter increments and row-address adds and
stores are counted.

In the affine table for this direction set
F1=alpha*e^2, F2=alpha*e and F3=alpha*e^3,
using a square, one cube and constant maps.
For each literal coefficient l initialize
LIN=theta^l*F1 and QUAD=theta^(2*l)*F2 by
APPLY or COPY. Emit four literal slices for
that coefficient with widths 13,13,13,11.
No slice has padding bits. Each uses 2^width
entries of ten words, starting at the running
slice pointer. Initialize its entry zero as
zero, except that the first slice of a
direction starts with F3. Advance the slice
pointer by 10*2^width after each slice.

For each literal bit r in the slice, read
the ten words of LIN and QUAD and XOR the
pairs to retain its column in registers.
The column before the following shifts is
its contribution to the affine map.
Starting k at zero, emit a runtime loop up
to k=2^r-1 with comparison and branch, increment
and jump. For each k compute addresses of
entry k and entry k+2^r by mul and add using
the slice base and a literal destination offset.
For each of the ten words, add the index to
the source address, load, XOR the column
word, add the index to the destination
address and store. This expands all subset
XOR rows by induction over r. Between bit
actions advance LIN by X and QUAD by X^2
with componentwise BSHIFT_1, including loads
and stores. Their initial theta factors
give the basis direction beta's column
alpha*(e*beta^2+e^2*beta). Together with
the first slice's constant, the table
therefore equals affine_entry for all keys.

Each copied entry costs at most 100 including
runtime loop control and its word operations.
There are at most 8191 copied entries per slice.
The literal column setup, shifts, initial row
and slice maintenance cost at most 5000, so
the per-slice bound is 825000. Outside the
40 slices for each direction, the five BMUL,
two field products, square and constant maps,
including at most eighteen theta maps,
fit under 140000 operations. Thus the counted
outer work is bounded by
500*(40*825000+140000)+1500000,
below 20000000000, including all explicit
branches and early exits.

At the inner VISIT the four V fields are
already populated. Guard and score cost
at most 425. On descent before changing
KEY1, save u,c,d at SAVED+30*depth using
the new depth. Emit one mul, one add and
one add, load and store per literal word,
cost 92. Apply the current direction's
affine slot to the old u, accumulating
into the ten old d words. Use mul and add for the direction address
and load the old ten d words. Per slice load
its u coefficient, shift and mask to the
literal width, multiply by ten and add the
literal offset and the direction base.
Per word add the index, load and XOR into its
accumulator. Store the ten updated words.
The cost is 2+20+40*(6+30)=1462. Then XOR the
two direction fields e,v*e into u,c from
DIRECTIONS+20*i, using mul and add and
add, two loads, XOR and store per word
for 102 primitives. FLIP changes KEY1.
On BACK load u,c,d from the current save
row before decrementing depth, via the same
92-operation row setup and transfer, then
undo FLIP. The maintained d is correct by
the cubic difference identity, and u,c
are correct by their linear relations.
The total inner visit plus entry and
exit updates costs at most
425+92+1462+102+92=2173, below 2200.
The ordinary key traversal is counted
separately. The field scratch and the
retained chart values use disjoint named
ranges. At most nineteen field accumulators and
fewer than thirty reused arithmetic or
address registers are live in any macro.
Four traversal registers and the chart
descriptor pointer, chart counter and
internal table counters fit in 64 registers
with those bounds.

## Nested key traversal and final output

Extract each of the forty received coefficients
into its literal Y array slot. Load the word
containing its first bit and shift right by its
literal offset. On crossing a word boundary,
load the next word, shift left and OR. Mask to
50 bits and store. This uses at most seven
operations each and at most 280 total. Set
BestD to 2001 by a move.

For the six ordinary charts emit literal
versions in pair order. Then emit one split
chart loop for each orientation, with the
v,u array identities fixed in that loop. In
each loop the chart counter runs from zero
to 63503, with explicit comparison, branch,
increment and jump. The descriptor pointer
starts at the orientation's first row and
increases by 16 after each chart. Its offsets
used during outer setup are literal. These
loops use one split chart body per orientation;
the Cartesian choices are stored as descriptor
data. A chart starts by copying its received
outer field to KEY0, setting depth to zero
and storing zero to NEXT[0]. This costs 21,
plus the split loop bookkeeping, all well
within 200 per split chart.

Each outer visit does the indicated chart
body. If its guard fails or the split matrix
is singular, branch to the outer EXTEND.
Otherwise an ordinary outer visit copies its
received inner group to KEY1, while split
outer setup has already placed the received
coefficients there. Initialize inner depth
to zero, compute its cap as chart cap minus
outer depth, and store zero to INEXT[0]. The
inner start outside precomputation costs at
most 22 including the possible ordinary copy.

Both traversals use the following template,
using their own KEY, depth, cap, NEXT, LAST
and labels. The inner DONE transfers to
outer EXTEND. The outer DONE falls through
to the next chart or split chart loop update.

    VISIT:
        execute the corresponding inlined visit
    EXTEND:
        flag := eq(depth,cap)
        branch if flag to BACK
        adr := add(NEXT,depth)
        i := load(adr)
        flag := eq(i,500)
        branch if flag to BACK
        j := add(i,1)
        store(adr,j)
        depth := add(depth,1)
        adr := add(NEXT,depth)
        store(adr,j)
        adr := add(LAST,depth)
        store(adr,i)
        [in a split inner traversal, save and update]
        FLIP(i)
        jump to VISIT
    BACK:
        flag := eq(depth,0)
        branch if flag to DONE
        adr := add(LAST,depth)
        i := load(adr)
        [in a split inner traversal, restore]
        depth := sub(depth,1)
        FLIP(i)
        jump to EXTEND

A branch not taken falls through. FLIP is

    word := div(i,50)
    bit := mod(i,50)
    mask := shl(1,bit)
    adr := add(KEY,word)
    v := load(adr)
    v := xor(v,mask)
    store(adr,v)

The two divisions use unsigned word instructions.
The indices stay in 0..499, so the flip preserves
the padding convention. A successful extension
costs 21, a failed one at most six, a nonroot
BACK 13 and a root BACK two, excluding the split
updates charged to the inner visit bound.
For any complete traversal these sum to at most
40 per node. The parent advances its NEXT and
initializes the child's NEXT to the next greater
bit, so a tree of cap r enumerates precisely all
subsets of its 500 positions of size at most r.
The LAST entry identifies the flip to undo as
it returns. Each chart has at most
O_r=sum_(i<=r) binom(500,i) outer visits.
Across all outer keys the total potential inner
visits is N_r=sum_(i<=r) binom(1000,i), by the
partition into two disjoint keys and their combined
weight bound. Failed outer guards only reduce this
count.

Pack the ten coefficient pieces of BEST0 followed
by the ten pieces of BEST1 into the output bits in
increasing order. For each piece load it. If aligned
store directly to the output word. Otherwise shift
left by the offset and load, OR and store there.
On crossing a word boundary also shift the original
piece right by 64 minus the offset and load, OR and
store in the next word. This costs at most nine
per piece, at most 180 overall. Assigned bit intervals
are disjoint. A direct aligned store only reaches an
output word without previously assigned bits, so
input order preserves all earlier pieces. Every
output word is written and the high 24 padding bits
of the last word are zero.

Fixed initialization, ordinary chart starts and final
output fit within 2000 operations, including final
loop exits. Per chart add (C+40)*N_r+20000000062*O_r,
using its stated inner bound C. The second term includes
outer setup, 22 inner initialization and 40 outer
traversal overhead. For each of the 127008 split charts,
count another 200 for its setup and loop bookkeeping.
The exact integer total is W in the manifest. The
declared computation checks W^125<=2^61228, the
memory extent, the tower and minor witnesses, exact
noise and error counts, and representative chart and
cube-update checks. The proof bounds the complete
search and all code paths analytically.
