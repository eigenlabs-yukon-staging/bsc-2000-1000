# Bilinear parities with six million split information charts

## Encoder and decoder

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
order in the F polynomial basis.

Let Q=q^2. The square root in K is
sqrt(x)=x^(2^499), a binary linear map. From the two
systematic field coordinates a,b the encoder sends

    (a,b,c,d),   c=a*b,   d=sqrt(a)*b.

There are exactly 2^1000 distinct codewords, with
systematic coordinates carried in the first 1000 bits.
Each group uses ten consecutive 50-bit coefficients.

The decoder includes a fallback candidate formed by
copying the two received systematic fields, replacing
a zero field by 1. It also searches the following
split charts. Choose either full known systematic
field v, four coefficients of the other systematic
field u, and three coefficients in each parity c,d.
Together these are 1000 known bits, with combined
distance cap 86. For a fixed v, solve the selected
chart equations if the recovery matrix is invertible.
Retain nonzero a,b only. Score every obtained candidate
and the fallback by full received Hamming distance, with
ties to smaller a+2^500*b in the specified bit order.
This gives a unique output on every input. The number
of choices per orientation is
binom(10,4)*binom(10,3)^2=3024000, hence 6048000 charts.

## Split chart equations

For known v=a, both c and d are F-linear maps of u=b.
Their matrices have columns v*theta^i and sqrt(v)*theta^i
in the selected basis. For known v=b, u=a, c has the
same multiplication columns v*theta^i. The d columns
are coefficients of v*sqrt(theta^i), applied to the
entrywise square roots of the u coefficients. For each
selected d row, square the received base coefficient
and the corresponding matrix row entries. This gives
an equivalent F-linear equation in the u coefficients,
because squaring is additive and bijective in F.

With the known four coefficients P and unknown six U,
the selected parity equations give M*U+J*P=Z.
Here the three equations from each parity are kept in
increasing row-index order, and Z is the observed six
parity coefficients with d transformed as stated when
needed. If M is invertible, U=M^-1*J*P+M^-1*Z.
Thus every eligible message can be recovered from its
known key, and each considered nonzero result exactly
matches the enumerated key. All matrix tests and the
remainder of the counted decoder are specified below.

## Regular messages and witness minors

In the a-outer orientation parameterize v by x=sqrt(v),
which is a bijection on K. Each of the three c matrix
rows is linear in x^2, hence has degree at most two in
the ten F coefficients of x. The three d rows are
linear in x. For b outer, parameterize by v directly.
The three c rows are linear and the three squared d
rows are at most quadratic in its coefficients. Thus
each 6-by-6 submatrix determinant polynomial has total
degree at most nine. Both declared minors runs verify
every determinant value is nonzero at the literal
witness for that orientation. This proves every
determinant polynomial is itself nonzero.

The witness calculation forms the twenty rows from
the two parity matrices directly. For fixed row and
column index sets of size k, expand the determinant
along the largest row, summing that row's entry times
the preceding row set's minors with one column
removed. Signs are plus in characteristic two. The
minors.py calculation uses size-k column masks,
looping over each c triple, all d prefixes of length
one and two, then all d triples. It computes all
210 minors of size six for each resulting row set.
Multiplication by a literal witness entry is computed
by XOR of seven byte-lookup contributions in F.
These tables are derived from its first fifty binary
basis multiples. This is an exact evaluation procedure
and both runs assert precisely 3024000 nonzero minors
and zero vanishing ones.

Over a field of size U, a nonzero polynomial
in n variables of total degree at most s has at
most s*U^(n-1) zeros. For n=1 this is the
univariate root bound. For larger n, write the
polynomial in its last variable with degree k
and leading coefficient of degree at most s-k
in the other variables. Count at most k choices
per prefix where the coefficient is nonzero and
use induction to count prefixes where it vanishes,
with at most U choices for the last variable there.
This gives the stated bound. Define a message as
regular if both fields meet every chart matrix
test in their own orientation. A zero field fails
the tests. The union bound over the 3024000 polynomials
gives at least q-27216000*2^450 possible values in each
orientation. The number of regular messages is
at least G=(q-27216000*2^450)^2.

## Conditional error bounds

Fix a noise vector z and its total weight t. For a
regular message, if its noise meets any split cap,
the true message is considered by the decoder.
It is then correct unless a distinct candidate
has full distance at most t. Counting all such
competitors from the code only enlarges this event.
We use three disjoint direction types according
to whether either systematic field is unchanged.
The bounds average over all uniform first messages
unless otherwise indicated, before restricting
to the regular subset.

For directions (s,u) with both nonzero, analyze
the second parity by squaring it first, giving
a*b^2. This is an invertible F2-linear change
of parity coordinates. The two differential
forms as functions of the first message are

    a*u+b*s+s*u,
    a*u^2+s*b^2+s*u^2.

They are affine binary forms. Fix a value of the
first. Solving for b leaves the second affine
in a with linear part (u^2/s)*(a^2+s*a).
Its kernel has exactly two points, 0 and s.
The joint distribution is therefore uniform on
an affine hyperplane in the 1000 binary parity
bits, with mass 2/Q on each point. The linear
bijective square-root map back to the actual d
preserves this property.

There is at most one such double direction per
nonzero binary hyperplane mask. The binary trace
Tr(x)=sum_(i=0..499) x^(2^i) is a nonzero F2-linear
function into F2 by the root bound. Multiplication
by a nonzero field scalar is surjective, so the
trace pairing is nonsingular and every binary
linear form on K is Tr(g*x) for a unique g.
A parity mask in the squared representation has
form Tr(alpha*c+beta*d_squared). To vanish on
the linear differential it must satisfy

    alpha*u+beta*u^2=0,
    alpha*s+sqrt(beta*s)=0.

The second uses Tr(x*b^2)=Tr(sqrt(x)*b).
With s,u nonzero and a nonzero mask, beta and
alpha are both nonzero. The equations give
u=alpha/beta and s=beta/alpha^2. This proves
the multiplicity bound, also in the actual
binary parity representation by invertibility.

For a fixed z and a competing systematic noise
of weight j, its direction is fixed by XOR.
The allowed parity noise lies in a 1000-bit
ball of radius ell=t-j. Denote its volume by J.
For a binary mask with at least D0=42 ones and
zeros, the signed shell count of weight k is the
coefficient of z^k in (1-z)^w(1+z)^(1000-w),
where w is the mask weight and the signs come
from the binary character. Pair D0 factors of
opposite signs and take absolute values of
coefficients. This bounds the absolute signed
ball count by H, the sum through ell of the
coefficients of (1+z^2)^D0*(1+z)^(1000-2*D0).
Translation by fixed actual parity noise changes
the signed sum by at most a sign. The affine
hyperplane therefore meets the translated ball
in at most (J+H)/2 points. H<=J by coefficient
comparison. If neither 500-bit half of the mask
has at least D0 ones and zeros, use the trivial
intersection bound J. The number of these
directions is at most Kbad=(2*sum_(i<D0) binom(500,i))^2
by the injective mask count. Thus among the
N=binom(1000,j) possible competing systematic
patterns at most min(N,Kbad) require the trivial
bound. Summing N*(J+H)+min(N,Kbad)*(J-H) over
j=0..t defines V_double(t). The expected count
of double-direction competitors is at most
V_double(t)/Q. Including patterns that do not
change both system fields only enlarges the
upper estimate.

Now suppose a is unchanged and b changes by
u!=0, with a nonzero as on the regular subset.
The three differences in the affected groups
are u, a*u, sqrt(a)*u. They are all nonzero and
any two determine the third and the unchanged
field a. For instance the ratio of the two
parity differences is sqrt(a). Similarly if b
is unchanged and a changes by s!=0, the three
differences are s, s*b, sqrt(s)*b. The ratio
of the parity differences is sqrt(s) and again
any two fix the whole triple and b. For fixed
z, suppose the competing noise in those three
groups has weights i,j,k. The numbers of
possible differences in their translated sets
are binom(500,i), binom(500,j), binom(500,k).
The number of valid triples is at most the
product of the two smallest of these counts,
since each projection onto two components
is injective. Each triple corresponds to the
unique unchanged field and the given nonzero
change. With the unchanged field uniform, the
expected competitor count contribution is at
most that product divided by q. We need this
count only on regular first messages, so we
may ignore the zero unchanged field and count
over all other possibilities.

Let V_partial(w) sum the products of two
smallest counts over i+j+k<=w. If the unchanged
group's actual noise has weight A, the competitor
has the same noise there, so use w=t-A. The
other orientation uses w=t-B. The computation
only needs w through 168, less than 250, so
binom(500,i) is increasing throughout these
shells. It sums sorted triples i<=j<=k with
the appropriate multiplicity, multiplying
binom(500,i)*binom(500,j). It confirms the
first crossing of q at w=169. This upper count
is nondecreasing in w.

Consequently for a noise z satisfying some
split cap, with systematic group weights
A,B, the success probability over uniform
messages is bounded below by

    max(0, G - V_double(t)
             - q*(V_partial(t-A)+V_partial(t-B))) / Q.

Indeed at least G regular messages are
eligible, and on each the decoder succeeds
unless there is a competitor of distance at
most t. The union bound on that event is
at most the expected count of competitors,
with the two partial terms estimated over
nonzero unchanged fields as above. This
argument covers the fallback as well as
every split candidate. When V_double(t)>=Q
or either partial bound reaches q we may
take zero, and these counts are nondecreasing
at the respective weights. The exact
certificate uses the first global crossing
after t=221 and the one-side crossings after
w=168. Noise vectors outside the counted
region contribute nonnegative success.

## Exact admitted noise count

For regular messages admission depends solely on the
error vector. Let A,B be its full systematic weights.
Let v4(X) be the sum of four smallest chunk weights of
a whole group X, and v3(X) the analogous sum of three
smallest. Then the minimum total chart noise weight is

    min(A+v4(B),B+v4(A)) + v3(C)+v3(D).

The symbols B or A inside v4 refer to the group, not
its total. All chunk choices are present and the
regular-message tests pass. Thus this minimum being
at most 86 specifies the exact needed condition.

For a single 500-bit group the checker obtains the
joint distribution by total weight w and sum of the
k lightest chunks, for k=3 and 4, through w<=221.
Let r be the kth lightest chunk weight. There are j<k
chunks strictly below r, e with weight r, and
10-j-e strictly above r, where e>=k-j. Choose their
positions by binom(10,j)*binom(10-j,e). The low and
high chunk total weights are enumerated by the
convolution powers of the polynomials with coefficients
binom(50,i) at weights i<r and i>r respectively.
The equal weights contribute binom(50,r)^e. If the low
and high totals are L,H, respectively, then w=L+e*r+H
and the k-lightest sum is L+(k-j)*r. Since at least
11-k chunks have weight at least r, it suffices to
iterate r<=floor(221/(11-k)). This partitions all
chunk weight outcomes exactly. The checker verifies
the full histogram mass binom(500,w).

Take the suffix in the k=4 histogram to obtain
J_w(r), the number of masks of total w whose four
lightest chunk weights sum to at least r. For a
fixed A<=B, the number with
u=min(A+v4(B),B+v4(A))>=j is
J_A(j-B)*J_B(j-A). Differences of neighboring
thresholds give equality counts. The threshold
range from A through
min(A+floor(4*B/10),B+floor(4*A/10)) contains every
u because the smallest four weights sum to at most
four tenths of the total. Double counts for A<B.

For total parity weight p, convolve the two k=3
histograms over the ordered totals C,D with C+D=p.
This gives the distribution of v3(C)+v3(D).
The verifier packs polynomial coefficients into
1024-bit slots to sum these convolutions as
exact integers. Every coefficient is at most
2^1000, so no cross-slot carries occur. It extracts
coefficients and checks the total binom(1000,p),
then takes their prefix sums. For each A,B,p with
t=A+B+p<=221 the number of admitted errors uses
the systematic equality classes and that parity
prefix at 86-u. Multiply the exact admitted
count by the nonnegative conditional success
numerator derived above. Sum those products
weighted by 9^(2000-t) and divide by Q*10^2000.
This lower-bounds average success. The exact
arithmetic checks that the complementary error
is at most 20939949382624787/125000000000000000.

## Scalar RAM representation and arithmetic

Use ten consecutive words for each K element, low
coefficient first, one 50-bit coefficient per word.
The unused high bits of each word are zero. The
62 ten-word arrays in ARRAYS start at 48. NEXT,
LAST,INEXT,ILAST have 87 words each, starting at
668,755,842,929. MULT has sixteen single-word rows
from 1016. The two reduction tables have two and
sixteen words at 1032 and 1034. GAMMA_MAP at 1050
contains two maps of seven slices of 256 one-word
entries, for multiplication by 184 and 108, as
given by gamma_entry. MAP_BASE at 4634 has the
19 MAP_TAGS maps, each 70 byte slices with 256
entries of ten words per slice. They implement
square root and the constant products theta^j
and sqrt(theta^j) for j=1..9; the literal j=0
map is COPY. table_columns and table_words fix
every entry including zero contributions from
unused high bits.

CHART_META at 3409434 stores the 6048000 descriptor
rows of 17 words each. descriptor(ordinal) determines
the table and lists orientation, four known indices,
six missing indices, and the two three-row sets.
The known-four and parity-three combinations are
in increasing lexicographic order. No host
experiment needs to materialize this entire
read-only table. DIRECTIONS at 106225434 has
500 writable rows of 30 words. All memory
ranges are disjoint and the largest exclusive
address is 106240434, below 2^32.

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
table_columns and table_words give the entire table.
The count is 70*(5+30)+10=2460. Each instruction is scalar.
Input and output may alias after the final read.
The same method for either map in GAMMA_MAP
costs at most 7*(5+3)+1=57 for a base element.

For MUL over K use schoolbook polynomial
convolution followed by descending modular
reduction. Use 19 literal coefficient slots in
PROD0 and PROD1. Store zero in every slot.
For each of the 100 ordered coefficient pairs
use BMUL into PRODTEMP, then XOR its result
into the appropriate slot (two loads, XOR and
store). For k=18 down through 10, replace
the contribution in position k by its terms
in positions k-7, k-9 and k-10, multiplying
by 1,108,184 respectively via GAMMA_MAP and
performing the three accumulator XOR updates.
Those positions are all lower, including when
they still require reduction on a later step.
Copy the lowest ten coefficients to the
destination. All field operands are distinct
from the product scratch, with aliasing of the
output allowed after the final input read.
The cost is at most
19+100*(251+4)+9*(2*57+3*4)+20=26673,
so use 30000 as a bound.

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

## Split chart precomputation and updates

For each orientation the code has one literal split
chart loop, with v,u array choices fixed. The runtime
descriptor pointer selects seventeen fixed-offset
indices for the current chart. The known and missing
coefficient indices are among 0..9.

At each outer visit copy KEY0 to v and guard zero.
Build the c multiplication columns into MC0..MC9
by COPY and APPLY of the theta powers. In a outer,
use APPLY of square root to put sqrt(v) into VROOT
and form MD0..MD9 from it by the same powers.
In b outer, use v and APPLY of the square roots
of the theta powers to build those columns. Each
relevant MD coefficient is then squared by BMUL
when taking it as an entry of M or J, only in b
outer. A column address is formed by multiplying
its descriptor index by ten and adding the fixed
column-array base. A row index is then added to
select one coefficient.

At GAUSS0 begin a 72-word augmented workspace for
[M|I], with six twelve-word rows. It fits within the
eight consecutive GAUSS arrays. Read M's six columns
and rows by literal descriptor offsets and store
them in the left half; use literal zero/one stores
for the identity. For each literal pivot column
choose the first nonzero entry in rows at or below
the diagonal, using load, eq and branch. If no
pivot exists skip this outer key. Otherwise swap
that row into the current row if needed, by load
and store pairs at each of twelve offsets. Use
the saved pivot in ROW_TMP to compute its inverse
into BASEINV by the binary exponentiation chain
for 2^50-2, with at most 100 BMUL. Multiply each
entry of the pivot row by the inverse using BMUL.
For every other literal row, preserve its pivot
entry and update each of its twelve entries by
XOR with that entry times the corresponding
normalized pivot row entry. The pivot process
reaches [I|M^-1] exactly when M is invertible.
All the addressed movement and comparisons in
this six-column loop are counted.

Write the six columns of M^-1 into COEF4..COEF9,
first six words of each, and the four columns
of M^-1*J into COEF0..COEF3 by six base products
and XOR for each scalar entry. The row ordering is
the missing coefficient order in the descriptor.
From the Y arrays read the observed four u chunks,
three c and three d chunks. Store the genuine
observations in KEY1, and the transformed values
in KEY_TRANSFORM, with the three d coefficients
squared only for b outer. Populate u at the root
by COPY of known values and dot products of all
ten transformed observations with the appropriate
coefficient entries for each missing value.
Compute c by MUL, and d by MUL of b and the square
root of a by APPLY. Thus V0..V3 initially
hold the root of the current inner search.

To build the direction for a key bit i=50*lane+bit,
start with f=X^bit, a single word obtained by
shl(1,bit). Emit ten literal lane bodies, each
with a fixed 50-iteration bit counter. For b outer
and lanes 7..9 use BMUL(f,f) to transform the
right-hand-side change. For each missing u coefficient
use BMUL by the corresponding COEF[lane] entry.
Include f in the changed known u coefficient if
lane<4 and zeros in the other known positions.
The full direction e goes into DELTA. Compute
DC=MUL(v,e). For a outer use DD=MUL(VROOT,e); for
b outer use APPLY(sqrt,e,TEMP) then DD=MUL(v,TEMP).
Each row at DIRECTIONS+30*i stores e,DC,DD in
order. Scalar address formation, copies and
branches for loop tests, increments and jumps
are included.

The initial matrix, inverse, coefficient and
root calculations use at most 1500000 operations.
For instance the six inverses use at most 600 BMUL,
normalization and elimination at most 432 BMUL,
and the coefficient and root dot products at
most 204 more. Other work includes the field maps,
a few full field products and literal word movements.
Each direction uses at most seven BMUL, two MUL,
one APPLY, and less than 3000 other operations
including index arithmetic, loop control and row
stores. The direction bound is below 70000,
so 500 of these and the earlier operations
cost below 50000000 per outer visit.

Each inner visit has a complete four-field candidate.
Use the two nonzero filters and SCORE at most 425.
On a descent or matching ascent, XOR the direction
fields into u,c,d. To do this emit mul and add for
the row base and per word add, two loads, XOR
and store. This costs 2+30*5=152 each time. The
solve is binary linear in key bits, including
the base squaring in b-outer d lanes, and both
parity functions are binary linear in u with
v fixed. Thus the direction is exactly the change
for the current key flip and its own inverse.
The ordinary key FLIP follows that update. The
invariant starts at the root and is maintained
under both descent and ascent. Including both
updates and a visit costs at most 425+2*152=729,
below 750 per inner node, in addition to the
traversal overhead charged below.

The augmented workspace, field product scratch,
coefficient arrays and retained chart arrays are
all disjoint in the listed layout. Within a macro
at most ten field accumulators plus twenty
word arithmetic and address registers are
needed simultaneously. BestD, the two depths,
inner cap, loop counters and descriptor pointer
keep the total below 64.


## Key traversal and output

Extract the forty received 50-bit chunks into Y0..Y3.
Each literal extraction loads its first input word and
shifts right by the fixed offset. When crossing a word,
load the next, shift left and OR. Mask to fifty bits and
store to the field array. This costs at most 280 in all.
Set BestD=2001, copy Y0,Y1 into V0,V1 and guard each zero,
replacing it with a literal 1 in its low word on zero.
Compute c by MUL and d by an APPLY(sqrt) and MUL, then
score the fallback. This fixed work, including both
field products, is below 70000.

For each of the two orientation loops, initialize its
descriptor pointer to that orientation's first row.
Visit its 3024000 rows in order, advancing the pointer
by seventeen and incrementing a counter. Use explicit
comparison, branch, increment and jump for the loop.
For each chart copy its received outer field to KEY0,
set outer depth zero and store zero at NEXT[0].
These instructions are covered by 200 per chart.
For each outer node execute the precomputation above,
including KEY1 at the inner root. On success set inner
depth zero, compute cap as 86 minus outer depth and
store zero at INEXT[0]. This costs two primitives, below
the allowance of 22 outside the outer body. A failed
outer guard skips the inner traversal.

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
        [in a split inner traversal, XOR the direction into u,c,d]
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

The base work outside the chart loops, including their
final comparisons, is below 100000. Per chart the upper
bound is (750+40)*N_r+50000062*O_r+200, where r=86,
and the outer term includes the setup, inner start and
outer traversal overhead. Multiplying by 6048000 and
adding the base work gives the stated integer W. The
certificate checks W^1000<=2^450797, both field
irreducibility certificates, exact error counts and
sample key reconstructions and direction updates.
Two further declared runs verify all rank witnesses.
All input and control paths are covered by the
bounded analytic program, without executing the full
candidate enumeration.
