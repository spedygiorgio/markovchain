# Restrict a Markov chain to a subset of states

Restricts a `markovchain` object to a chosen subset of its states,
either as a raw (generally non-stochastic) principal submatrix, or as a
properly renormalized Markov chain describing behaviour conditional on
staying inside the subset.

## Usage

``` r
subchain(object, states, method = c("submatrix", "renormalize"))

# S4 method for class 'markovchain'
subchain(object, states, method = c("submatrix", "renormalize"))
```

## Arguments

- object:

  A `markovchain` object.

- states:

  A character vector of at least one state name from `states(object)`,
  with no duplicates: the subset to restrict to.

- method:

  Either `"submatrix"` or `"renormalize"` (can be abbreviated). See
  Details: the two options answer genuinely different questions and are
  not interchangeable.

## Value

If `method = "submatrix"`: a plain numeric matrix (*not* a `markovchain`
object, since its rows generally do not sum to one), the principal
submatrix of the transition matrix restricted to `states`, always
returned in row-stochastic orientation regardless of `object`'s own
storage convention.

If `method = "renormalize"`: a new `markovchain` object on exactly the
states in `states`, row-stochastic, describing the chain conditional on
never leaving that subset.

## Details

The two methods are deliberately named after two different, standard
constructions, so that the choice – and its consequences – is explicit
rather than implied:

- `"submatrix"`:

  Simply the entries of \\P\\ with both indices restricted to `states`,
  with no adjustment. Its rows generally sum to *less than* one, because
  probability mass that originally went to states outside the subset is
  dropped, not redistributed. This is the "\\Q\\" block used, e.g., when
  building the fundamental matrix of an absorbing chain (see
  [`fundamentalMatrix`](fundamentalMatrix.md)): a useful building block
  for other computations, but not itself a transition matrix of any
  Markov chain, which is why it is returned as a plain matrix.

- `"renormalize"`:

  Each retained row is divided by its own sum, so the result is
  row-stochastic and can be wrapped in a `markovchain` object. This is
  the chain of successive positions of `object`, *conditioned on the
  event that it never leaves `states`* (sometimes called the chain
  "watched on" `states`, or its taboo probabilities; see Norris (1997),
  Section 3.3). It requires every state in `states` to have strictly
  positive probability of transitioning within the subset (otherwise
  that conditioning event has probability zero from that state, and the
  row cannot be renormalized); an error is raised naming any state that
  fails this, rather than silently producing a row of `NaN`.

Neither method requires `object` to be irreducible: restricting to a
subset of states is meaningful for any chain, and is often used
precisely to study one communicating class in isolation.

The implementation performs no eigendecomposition; it is \\O(k^2)\\ time
and memory for a subset of size \\k\\, after an \\O(n^2)\\ extraction
from the full \\n\\-state matrix.

## References

Norris, J. R. (1998). *Markov Chains*. Cambridge University Press.

## See also

[`lazyChain`](lazyChain.md),
[`fundamentalMatrix`](fundamentalMatrix.md),
[`canonicForm`](structuralAnalysis.md)

## Examples

``` r
statesNames <- c("a", "b", "c")
mc <- new("markovchain", states = statesNames,
  transitionMatrix = matrix(c(0.5, 0.3, 0.2,
                              0.2, 0.6, 0.2,
                              0.1, 0.1, 0.8), byrow = TRUE, nrow = 3,
                            dimnames = list(statesNames, statesNames)))

# Raw submatrix: rows no longer sum to 1, mass has "leaked" to "c".
subchain(mc, c("a", "b"), method = "submatrix")
#>     a   b
#> a 0.5 0.3
#> b 0.2 0.6
rowSums(subchain(mc, c("a", "b"), method = "submatrix"))
#>   a   b 
#> 0.8 0.8 

# Renormalized: a genuine markovchain, conditional on staying in {a, b}.
watched <- subchain(mc, c("a", "b"), method = "renormalize")
watched
#> Unnamed Markov chain (subchain) 
#>  A  2 - dimensional discrete Markov Chain defined by the following states: 
#>  a, b 
#>  The transition matrix  (by rows)  is defined as follows: 
#>       a     b
#> a 0.625 0.375
#> b 0.250 0.750
#> 
rowSums(watched@transitionMatrix)
#> a b 
#> 1 1 
```
