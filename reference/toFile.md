# Write or read a Markov chain to or from a file

Writes a `markovchain` object to a JSON, YAML or CSV file, or reads one
back, using the same representation as
[`toDictionary`](toDictionary.md).

## Usage

``` r
toFile(object, file, format = NULL)

# S4 method for class 'markovchain'
toFile(object, file, format = NULL)

fromFile(file, format = NULL)
```

## Arguments

- object:

  A `markovchain` object (for `toFile`).

- file:

  A single file path to write to or read from. If `format` is not
  supplied, it is inferred from the file extension (`.json`,
  `.yaml`/`.yml` or `.csv`).

- format:

  One of `"json"`, `"yaml"` or `"csv"`. The default, `NULL`, infers the
  format from `file`'s extension.

## Value

`toFile` returns `file`, invisibly. `fromFile` returns a `markovchain`
object.

## Details

The JSON and YAML formats store exactly what
[`toDictionary`](toDictionary.md) returns (name, state names, and the
transition probabilities nested by source state), and round-trip a chain
exactly, including full numeric precision (`toFile` writes both JSON and
YAML with 17 significant digits, enough to recover every double
exactly).

The CSV format only stores the transition matrix itself, as a table of
probabilities with the state names as both the header row and the first
column – there is no natural place in a CSV file for the chain's `name`,
so it is not preserved by `toFile(..., format = "csv")` and `fromFile`
always returns an unnamed chain for a `.csv` file. This is the same
limitation PyDTMC's own CSV format has.

Writing JSON requires the jsonlite package, and writing YAML requires
the yaml package; both are only in `Suggests`, and an informative error
is raised if the relevant package is not installed. Reading has the same
requirements for the format being read. CSV uses only base R and has no
extra package dependency.

## See also

[`toDictionary`](toDictionary.md), [`fromDictionary`](toDictionary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
statesNames <- c("a", "b")
mc <- new("markovchain", states = statesNames,
          transitionMatrix = matrix(c(0.7, 0.3, 0.4, 0.6), byrow = TRUE,
                                     nrow = 2, dimnames = list(statesNames, statesNames)))
toFile(mc, "chain.json")
identical(fromFile("chain.json")@transitionMatrix, mc@transitionMatrix)
} # }
```
