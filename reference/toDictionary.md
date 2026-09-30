# Represent a Markov chain as a plain R list

Converts a `markovchain` object to a plain, self-describing R list: the
same information [`toFile`](toFile.md) writes to disk, kept in memory.
`fromDictionary` reverses the conversion.

## Usage

``` r
toDictionary(object)

# S4 method for class 'markovchain'
toDictionary(object)

fromDictionary(d)
```

## Arguments

- object:

  A `markovchain` object.

- d:

  A list as returned by `toDictionary`. `d$transitionMatrix` may also be
  a plain `n x n` matrix (with or without dimnames matching `d$states`)
  for convenience when building a dictionary by hand rather than from an
  existing `markovchain` object. When it is a plain matrix, `d$byrow` is
  honored exactly as `new("markovchain", ...)` honors its own `byrow`
  argument: the matrix is stored as given, with `d$byrow` only
  documenting whether it is row- or column-stochastic (so a
  column-stochastic matrix round-trips by setting `d$byrow = FALSE`,
  with no transposition performed here). When `d$transitionMatrix` is a
  nested list (as `toDictionary` produces), the nesting itself is always
  keyed `[[from]][[to]]` – i.e. row-stochastic – so it is always
  reconstructed with `byrow = TRUE`, regardless of `d$byrow`.

## Value

A named list with four elements:

- `name`:

  The chain's `name`, as a single string (possibly empty).

- `states`:

  A character vector of state names, in order.

- `byrow`:

  Always `TRUE`: the list always stores the chain row-stochastically,
  regardless of `object`'s own storage convention, so that the
  representation is unambiguous without also having to interpret this
  flag.

- `transitionMatrix`:

  A named list of named lists: `transitionMatrix[[i]][[j]]` is the
  probability of moving from state `i` to state `j`. This is
  deliberately not a plain matrix, so that the structure serializes to
  JSON or YAML (via [`toFile`](toFile.md)) as a self-describing object
  keyed by state name, rather than a bare array whose meaning depends on
  remembering a row/column order.

`fromDictionary` returns a `markovchain` object.

## Details

Unlike PyDTMC's own `to_dictionary()`/`from_dictionary()`, which
represent a chain as a flat mapping from every `(from_state, to_state)`
pair to its probability – \\n^2\\ entries with no state grouping – this
nests the representation by source state, which is both more compact to
read and directly round-trips through R's own list-of-lists idiom
without any special tuple-key handling.

## See also

[`toFile`](toFile.md), [`fromFile`](toFile.md)

## Examples

``` r
statesNames <- c("a", "b")
mc <- new("markovchain", states = statesNames,
          transitionMatrix = matrix(c(0.7, 0.3, 0.4, 0.6), byrow = TRUE,
                                     nrow = 2, dimnames = list(statesNames, statesNames)))
d <- toDictionary(mc)
d$transitionMatrix$a$b # 0.3: probability of moving from "a" to "b"
#> [1] 0.3

identical(fromDictionary(d)@transitionMatrix, mc@transitionMatrix)
#> [1] TRUE
```
