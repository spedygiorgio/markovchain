# Merge two Markov chains by convex combination of their transition matrices

Builds a new `markovchain` object whose transition matrix is a convex
combination of the transition matrices of two existing chains defined on
the same state space.

## Usage

``` r
mergeWith(object, other, gamma = 0.5)

# S4 method for class 'markovchain,markovchain'
mergeWith(object, other, gamma = 0.5)
```

## Arguments

- object:

  A `markovchain` object.

- other:

  A second `markovchain` object, defined on the same set of state names
  as `object` (see Details for what "same" means here).

- gamma:

  A single number in \\\[0,1\]\\: the weight given to `other`.
  `gamma = 0` returns `object` unchanged (up to storage convention);
  `gamma = 1` returns `other` unchanged.

## Value

A new `markovchain` object, row-stochastic, on the common state names,
with transition matrix \\P = (1-\gamma)P_1+\gamma P_2\\.

## Details

For transition matrices \\P_1\\ (from `object`) and \\P_2\\ (from
`other`), and a blending factor \\\gamma\in\[0,1\]\\, the merged
transition matrix is \$\$P = (1-\gamma) P_1 + \gamma P_2.\$\$

**States are matched by name, not by position.** `object` and `other`
must have exactly the same set of state names (as sets – `other`'s
states may be in a different order, or `other` may use a different
storage convention (`byrow`), and both are handled correctly). Rows and
columns of `other`'s transition matrix are realigned to `object`'s state
order before combining, and both matrices are converted to
row-stochastic form first if needed, so that \\P_1\\ and \\P_2\\ are
always combined entry-for-entry between matching states rather than
between matching matrix positions.

This is a deliberate difference from the na\\ive version of this
operation, which combines two same-*size* transition matrices
positionally and would silently produce a meaningless result if the two
chains happened to list their states in a different order (or under a
different `byrow` convention) despite describing the same states.
Requiring identical state name sets, rather than merely identical size,
catches that mismatch as an error instead of propagating it.

Because \\P_1\\ and \\P_2\\ are both row-stochastic with non-negative
entries and \\\gamma\in\[0,1\]\\, \\P\\ is automatically row-stochastic
with non-negative entries: no renormalization is needed (this is the
same convexity argument used for [`lazyChain`](lazyChain.md), of which
`mergeWith(object, identity_chain, gamma)` is a special case when
`other` is an identity chain on the same states).

This function does not require `object` or `other` to be irreducible:
merging is meaningful for any two chains on the same state space,
including reducible ones. Note, however, that the merged chain's
stationary distribution (if any) is generally *not* a combination of
\\\pi_1\\ and \\\pi_2\\ in any simple way; it must be recomputed from
\\P\\ directly.

The implementation performs no eigendecomposition and is \\O(n^2)\\ time
and memory for two \\n\\-state chains, dominated by realigning `other`'s
matrix to `object`'s state order.

## See also

[`lazyChain`](lazyChain.md), [`subchain`](subchain.md)

## Examples

``` r
statesNames <- c("a", "b")
mc1 <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0.9, 0.1, 0.1, 0.9), byrow = TRUE, nrow = 2,
                            dimnames = list(statesNames, statesNames)))
# mc2 lists the same two states in the opposite order.
mc2 <- new("markovchain", states = rev(statesNames),
  transitionMatrix = matrix(c(0.5, 0.5, 0.5, 0.5), byrow = TRUE, nrow = 2,
                            dimnames = list(rev(statesNames), rev(statesNames))))
merged <- mergeWith(mc1, mc2, gamma = 0.5)
merged
#> Unnamed Markov chain + Unnamed Markov chain (merged, gamma = 0.5) 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by rows)  is defined as follows: 
#>     a   b
#> a 0.7 0.3
#> b 0.3 0.7
#> 
# States are matched by name: merged["a", "b"] combines mc1["a","b"] with
# mc2["a","b"], not with whatever happened to sit in the same matrix cell.
```
