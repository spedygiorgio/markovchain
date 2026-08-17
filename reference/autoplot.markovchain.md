# Plot a Markov chain with ggplot2

Creates a ggplot2 representation of a discrete-time Markov chain. Nodes
are states and directed edges represent transitions with positive
probability. Communicating classes are shown with different node fills.

## Usage

``` r
autoplot.markovchain(
  object,
  threshold = 0,
  show_probabilities = TRUE,
  digits = 2,
  node_size = 6,
  edge_width = 1,
  ...
)
```

## Arguments

- object:

  An object of class \`markovchain\`.

- threshold:

  Minimum transition probability to display. Defaults to 0.

- show_probabilities:

  Logical; whether to label edges with transition probabilities.
  Defaults to \`TRUE\`.

- digits:

  Number of digits used to format transition probabilities.

- node_size:

  Size of state nodes in the ggplot2 plot.

- edge_width:

  Minimum width multiplier for transition edges.

- ...:

  Currently unused, reserved for future extensions.

## Value

A ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  weather <- matrix(c(0.7, 0.2, 0.1,
                      0.3, 0.4, 0.3,
                      0.2, 0.45, 0.35),
                    nrow = 3, byrow = TRUE,
                    dimnames = list(c("sunny", "cloudy", "rain"),
                                    c("sunny", "cloudy", "rain")))
  mc <- new("markovchain", states = rownames(weather),
            transitionMatrix = weather, name = "Weather")
  ggplot2::autoplot(mc)
}
#> Scale for linewidth is already present.
#> Adding another scale for linewidth, which will replace the existing scale.
```
