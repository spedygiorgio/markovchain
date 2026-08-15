# Function to generate a sequence of states from homogeneous or non-homogeneous Markov chains.

Provided any `markovchain` or `markovchainList` objects, it returns a
sequence of states coming from the underlying stationary distribution.

## Usage

``` r
rmarkovchain(
  n,
  object,
  what = "data.frame",
  useRCpp = TRUE,
  parallel = FALSE,
  num.cores = NULL,
  ...
)
```

## Arguments

- n:

  Sample size

- object:

  Either a `markovchain` or a `markovchainList` object

- what:

  It specifies whether either a `data.frame` or a `matrix` (each rows
  represent a simulation) or a `list` is returned.

- useRCpp:

  Boolean. Should RCpp fast implementation being used? Default is yes.

- parallel:

  Boolean. Should parallel implementation being used? Default is yes.

- num.cores:

  Number of Cores to be used

- ...:

  additional parameters passed to the internal sampler

## Value

Character Vector, data.frame, list or matrix

## Details

When a homogeneous process is assumed (`markovchain` object) a sequence
is sampled of size n. When a non - homogeneous process is assumed, n
samples are taken but the process is assumed to last from the begin to
the end of the non-homogeneous markov process.

## Note

Check the type of input

## References

A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall
2010

## See also

[`markovchainFit`](markovchainFit.md),
[`markovchainSequence`](markovchainSequence.md)

## Author

Giorgio Spedicato

## Examples

``` r
# define the markovchain object
statesNames <- c("a", "b", "c")
mcB <- new("markovchain", states = statesNames, 
   transitionMatrix = matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), 
   nrow = 3, byrow = TRUE, dimnames = list(statesNames, statesNames)))

# show the sequence
outs <- rmarkovchain(n = 100, object = mcB, what = "list")


#define markovchainList object
statesNames <- c("a", "b", "c")
mcA <- new("markovchain", states = statesNames, transitionMatrix = 
   matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
   byrow = TRUE, dimnames = list(statesNames, statesNames)))
mcB <- new("markovchain", states = statesNames, transitionMatrix = 
   matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
   byrow = TRUE, dimnames = list(statesNames, statesNames)))
mcC <- new("markovchain", states = statesNames, transitionMatrix = 
   matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
   byrow = TRUE, dimnames = list(statesNames, statesNames)))
mclist <- new("markovchainList", markovchains = list(mcA, mcB, mcC)) 

# show the list of sequence
rmarkovchain(100, mclist, "list")
#> [[1]]
#> [1] "b" "c" "b"
#> 
#> [[2]]
#> [1] "b" "c" "b"
#> 
#> [[3]]
#> [1] "c" "b" "c"
#> 
#> [[4]]
#> [1] "c" "b" "c"
#> 
#> [[5]]
#> [1] "c" "a" "c"
#> 
#> [[6]]
#> [1] "b" "c" "b"
#> 
#> [[7]]
#> [1] "b" "c" "b"
#> 
#> [[8]]
#> [1] "a" "b" "c"
#> 
#> [[9]]
#> [1] "b" "c" "b"
#> 
#> [[10]]
#> [1] "a" "b" "c"
#> 
#> [[11]]
#> [1] "c" "b" "b"
#> 
#> [[12]]
#> [1] "b" "c" "b"
#> 
#> [[13]]
#> [1] "b" "b" "c"
#> 
#> [[14]]
#> [1] "b" "b" "b"
#> 
#> [[15]]
#> [1] "c" "b" "b"
#> 
#> [[16]]
#> [1] "a" "a" "b"
#> 
#> [[17]]
#> [1] "c" "b" "c"
#> 
#> [[18]]
#> [1] "b" "b" "b"
#> 
#> [[19]]
#> [1] "b" "c" "b"
#> 
#> [[20]]
#> [1] "b" "c" "b"
#> 
#> [[21]]
#> [1] "c" "b" "c"
#> 
#> [[22]]
#> [1] "b" "c" "b"
#> 
#> [[23]]
#> [1] "b" "c" "b"
#> 
#> [[24]]
#> [1] "c" "b" "c"
#> 
#> [[25]]
#> [1] "c" "b" "c"
#> 
#> [[26]]
#> [1] "a" "c" "b"
#> 
#> [[27]]
#> [1] "c" "b" "b"
#> 
#> [[28]]
#> [1] "c" "b" "b"
#> 
#> [[29]]
#> [1] "b" "c" "b"
#> 
#> [[30]]
#> [1] "c" "b" "c"
#> 
#> [[31]]
#> [1] "c" "b" "c"
#> 
#> [[32]]
#> [1] "b" "c" "b"
#> 
#> [[33]]
#> [1] "b" "c" "b"
#> 
#> [[34]]
#> [1] "b" "c" "b"
#> 
#> [[35]]
#> [1] "c" "b" "c"
#> 
#> [[36]]
#> [1] "c" "b" "c"
#> 
#> [[37]]
#> [1] "c" "b" "c"
#> 
#> [[38]]
#> [1] "b" "b" "b"
#> 
#> [[39]]
#> [1] "b" "c" "c"
#> 
#> [[40]]
#> [1] "b" "b" "c"
#> 
#> [[41]]
#> [1] "b" "c" "b"
#> 
#> [[42]]
#> [1] "a" "c" "b"
#> 
#> [[43]]
#> [1] "b" "c" "b"
#> 
#> [[44]]
#> [1] "c" "b" "c"
#> 
#> [[45]]
#> [1] "a" "b" "c"
#> 
#> [[46]]
#> [1] "b" "b" "b"
#> 
#> [[47]]
#> [1] "b" "c" "b"
#> 
#> [[48]]
#> [1] "b" "c" "b"
#> 
#> [[49]]
#> [1] "b" "c" "a"
#> 
#> [[50]]
#> [1] "b" "c" "b"
#> 
#> [[51]]
#> [1] "c" "b" "c"
#> 
#> [[52]]
#> [1] "c" "c" "b"
#> 
#> [[53]]
#> [1] "b" "c" "b"
#> 
#> [[54]]
#> [1] "c" "b" "c"
#> 
#> [[55]]
#> [1] "a" "b" "b"
#> 
#> [[56]]
#> [1] "c" "b" "c"
#> 
#> [[57]]
#> [1] "b" "b" "c"
#> 
#> [[58]]
#> [1] "b" "c" "b"
#> 
#> [[59]]
#> [1] "c" "a" "b"
#> 
#> [[60]]
#> [1] "a" "c" "b"
#> 
#> [[61]]
#> [1] "b" "c" "b"
#> 
#> [[62]]
#> [1] "b" "c" "b"
#> 
#> [[63]]
#> [1] "b" "b" "c"
#> 
#> [[64]]
#> [1] "c" "b" "c"
#> 
#> [[65]]
#> [1] "c" "b" "c"
#> 
#> [[66]]
#> [1] "b" "b" "b"
#> 
#> [[67]]
#> [1] "c" "b" "c"
#> 
#> [[68]]
#> [1] "b" "c" "b"
#> 
#> [[69]]
#> [1] "b" "c" "b"
#> 
#> [[70]]
#> [1] "c" "a" "a"
#> 
#> [[71]]
#> [1] "b" "c" "a"
#> 
#> [[72]]
#> [1] "c" "b" "b"
#> 
#> [[73]]
#> [1] "b" "b" "c"
#> 
#> [[74]]
#> [1] "b" "c" "b"
#> 
#> [[75]]
#> [1] "c" "b" "c"
#> 
#> [[76]]
#> [1] "c" "b" "c"
#> 
#> [[77]]
#> [1] "c" "a" "c"
#> 
#> [[78]]
#> [1] "c" "b" "c"
#> 
#> [[79]]
#> [1] "c" "b" "c"
#> 
#> [[80]]
#> [1] "c" "c" "b"
#> 
#> [[81]]
#> [1] "c" "c" "b"
#> 
#> [[82]]
#> [1] "c" "b" "c"
#> 
#> [[83]]
#> [1] "b" "c" "b"
#> 
#> [[84]]
#> [1] "a" "a" "a"
#> 
#> [[85]]
#> [1] "c" "a" "c"
#> 
#> [[86]]
#> [1] "b" "c" "b"
#> 
#> [[87]]
#> [1] "b" "c" "b"
#> 
#> [[88]]
#> [1] "c" "b" "c"
#> 
#> [[89]]
#> [1] "c" "b" "c"
#> 
#> [[90]]
#> [1] "b" "c" "b"
#> 
#> [[91]]
#> [1] "b" "c" "b"
#> 
#> [[92]]
#> [1] "b" "c" "b"
#> 
#> [[93]]
#> [1] "c" "b" "b"
#> 
#> [[94]]
#> [1] "b" "c" "b"
#> 
#> [[95]]
#> [1] "c" "b" "c"
#> 
#> [[96]]
#> [1] "b" "b" "c"
#> 
#> [[97]]
#> [1] "b" "c" "b"
#> 
#> [[98]]
#> [1] "b" "c" "a"
#> 
#> [[99]]
#> [1] "a" "b" "c"
#> 
#> [[100]]
#> [1] "b" "c" "b"
#> 
     
```
