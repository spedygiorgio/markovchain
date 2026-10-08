# Export the transition graph as Graphviz DOT or Mermaid text

`toDot` writes the transition graph of a Markov chain in the DOT
language of Graphviz; `toMermaid` writes it as a Mermaid flowchart. The
states are the nodes and every transition with positive probability is
an edge labelled with its probability.

## Usage

``` r
toDot(object, file = NULL, digits = 3, minProbability = 0, direction = "LR")

# S4 method for class 'markovchain'
toDot(object, file = NULL, digits = 3, minProbability = 0, direction = "LR")

toMermaid(
  object,
  file = NULL,
  digits = 3,
  minProbability = 0,
  direction = "LR"
)

# S4 method for class 'markovchain'
toMermaid(
  object,
  file = NULL,
  digits = 3,
  minProbability = 0,
  direction = "LR"
)
```

## Arguments

- object:

  A `markovchain` object.

- file:

  An optional file path. If given, the text is also written there
  (UTF-8) and returned invisibly.

- digits:

  The number of significant digits of the probabilities shown on the
  edges.

- minProbability:

  Transitions with probability not above this value are left out; the
  default, 0, keeps every possible transition.

- direction:

  The layout direction: `"LR"` (left to right, the default), `"TB"` (top
  to bottom), `"RL"` or `"BT"`.

## Value

The text, as a single character string (invisibly when `file` is given).

## Details

The edges follow the outgoing distribution of each state whatever the
storage of the transition matrix (`byrow`), so a chain and its
column-stochastic copy give the same text. State names are quoted and
escaped, so they may contain spaces, quotes or other punctuation. In the
Mermaid output the nodes get the identifiers `s1`, `s2`, ..., with the
state names as labels, because Mermaid identifiers cannot contain
arbitrary characters.

The DOT text can be rendered with Graphviz (for instance
`dot -Tpng chain.dot -o chain.png`) or with
[`DiagrammeR::grViz()`](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html);
the Mermaid text can be pasted into any Markdown renderer that supports
Mermaid diagrams, or rendered with
[`DiagrammeR::mermaid()`](https://rich-iannone.github.io/DiagrammeR/reference/mermaid.html).

## See also

[`toFile`](toFile.md),
[`plot,markovchain,missing-method`](markovchain-class.md)

## Examples

``` r
statesNames <- c("a", "b", "c")
mc <- new("markovchain", states = statesNames, name = "Example",
  transitionMatrix = matrix(c(0.5, 0.5, 0, 0.2, 0.3, 0.5, 0, 0, 1),
    byrow = TRUE, nrow = 3, dimnames = list(statesNames, statesNames)))
cat(toDot(mc), "\n")
#> digraph "Example" {
#>   rankdir=LR;
#>   node [shape=circle];
#>   "a";
#>   "b";
#>   "c";
#>   "a" -> "a" [label="0.5"];
#>   "a" -> "b" [label="0.5"];
#>   "b" -> "a" [label="0.2"];
#>   "b" -> "b" [label="0.3"];
#>   "b" -> "c" [label="0.5"];
#>   "c" -> "c" [label="1"];
#> } 
cat(toMermaid(mc, direction = "TB"), "\n")
#> flowchart TB
#>   s1(("a"))
#>   s2(("b"))
#>   s3(("c"))
#>   s1 -->|"0.5"| s1
#>   s1 -->|"0.5"| s2
#>   s2 -->|"0.2"| s1
#>   s2 -->|"0.3"| s2
#>   s2 -->|"0.5"| s3
#>   s3 -->|"1"| s3 
```
