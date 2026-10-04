# Build a birth-death Markov chain

Constructs a `markovchain` object for a birth-death process: a chain on
linearly ordered states \\1,2,\ldots,n\\ that, from any state, can only
move to itself or to an immediately adjacent state.

## Usage

``` r
birthDeath(p, q, states = NULL)
```

## Arguments

- p:

  A numeric vector of length \\n-1\\: `p[i]` is the "birth" probability
  of moving from state \\i\\ up to state \\i+1\\.

- q:

  A numeric vector of length \\n-1\\, the same length as `p`: `q[i]` is
  the "death" probability of moving from state \\i+1\\ down to state
  \\i\\.

- states:

  An optional character vector of \\n=\code{length(p)}+1\\ state names.
  Defaults to `as.character(1:n)`.

## Value

A new, row-stochastic `markovchain` object on \\n\\ states, with
transition matrix \$\$P\_{ii}=1-p_i-q\_{i-1},\quad P\_{i,i+1}=p_i,\quad
P\_{i,i-1}=q\_{i-1}\$\$ (boundary terms \\q_0\\ and \\p_n\\ are
understood to not exist, i.e. \\P\_{11}=1-p_1\\ and
\\P\_{nn}=1-q\_{n-1}\\).

## Details

**Why `p` and `q` have length \\n-1\\, not \\n\\.** Every birth-death
transition is a move between two adjacent states, and there are exactly
\\n-1\\ adjacent pairs among \\n\\ linearly ordered states:
`p[i]`/`q[i]` unambiguously describe the pair \\(i,i+1)\\. This
sidesteps a common source of confusion in this construction, namely what
to do with a stray "birth probability of the top state" or "death
probability of the bottom state" – quantities that do not correspond to
any actual transition, since there is no state \\n+1\\ to be born into
or state \\0\\ to die into. Some implementations accept two length-\\n\\
vectors and quietly renormalize every row so that any such leftover
probability mass is redistributed among the transitions that do exist;
`birthDeath()` instead makes the \\n-1\\ genuine transition
probabilities the only inputs, so there is no leftover mass to
(silently) dispose of in the first place.

Every row's diagonal entry is determined by the requirement that the row
sums to \\1\\, so `p` and `q` alone fully determine \\P\\: no separate
"staying" probability is accepted or needed. `p+q` is allowed to reach
\\1\\ for an interior state (no staying probability there), but each
element of `p` and `q` must itself lie in \\\[0,1\]\\ and `p[i]+q[i]`
for the shared index \\i\\ need not be checked against 1 the way it
would for a single state's own two probabilities, since `p[i]` leaves
state \\i\\ while `q[i]` leaves state \\i+1\\: the actual per-state
constraint, \\p_i+q\_{i-1}\le 1\\, is checked directly on the assembled
diagonal.

The two boundary states \\1\\ and \\n\\ are reflecting only in the weak
sense that no birth/death carries them outside \\\\1,\ldots,n\\\\ – they
still generally have a positive probability of staying put (\\1-p_1\\
and \\1-q\_{n-1}\\ respectively) rather than being forced to bounce
back, unlike [`gamblersRuin`](gamblersRuin.md)'s absorbing ends or
[`toBoundedChain`](toBoundedChain.md)'s explicit reflecting condition,
which can be applied afterwards to force deterministic bouncing or
absorption at the ends of any chain, including one built here.

## See also

[`gamblersRuin`](gamblersRuin.md),
[`toBoundedChain`](toBoundedChain.md), [`urnModel`](urnModel.md)

## Examples

``` r
# A simple 4-state birth-death chain with constant birth/death rates.
bd <- birthDeath(p = c(0.3, 0.4, 0.5), q = c(0.2, 0.3, 0.1))
bd
#> Birth-Death Chain 
#>  A  4 - dimensional discrete Markov Chain defined by the following states: 
#>  1, 2, 3, 4 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     1   2   3   4
#> 1 0.7 0.3 0.0 0.0
#> 2 0.2 0.4 0.4 0.0
#> 3 0.0 0.3 0.2 0.5
#> 4 0.0 0.0 0.1 0.9
#> 
rowSums(bd@transitionMatrix)
#> 1 2 3 4 
#> 1 1 1 1 
```
