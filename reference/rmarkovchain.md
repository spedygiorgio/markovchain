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
#> [1] "b" "c" "a"
#> 
#> [[2]]
#> [1] "c" "b" "b"
#> 
#> [[3]]
#> [1] "c" "b" "c"
#> 
#> [[4]]
#> [1] "b" "c" "b"
#> 
#> [[5]]
#> [1] "c" "b" "c"
#> 
#> [[6]]
#> [1] "c" "b" "b"
#> 
#> [[7]]
#> [1] "a" "b" "c"
#> 
#> [[8]]
#> [1] "c" "b" "c"
#> 
#> [[9]]
#> [1] "b" "b" "c"
#> 
#> [[10]]
#> [1] "b" "b" "c"
#> 
#> [[11]]
#> [1] "c" "b" "c"
#> 
#> [[12]]
#> [1] "a" "b" "c"
#> 
#> [[13]]
#> [1] "c" "b" "c"
#> 
#> [[14]]
#> [1] "b" "b" "b"
#> 
#> [[15]]
#> [1] "b" "c" "b"
#> 
#> [[16]]
#> [1] "c" "b" "b"
#> 
#> [[17]]
#> [1] "c" "a" "b"
#> 
#> [[18]]
#> [1] "a" "c" "b"
#> 
#> [[19]]
#> [1] "a" "b" "c"
#> 
#> [[20]]
#> [1] "b" "b" "b"
#> 
#> [[21]]
#> [1] "c" "b" "c"
#> 
#> [[22]]
#> [1] "b" "b" "c"
#> 
#> [[23]]
#> [1] "b" "b" "c"
#> 
#> [[24]]
#> [1] "c" "a" "c"
#> 
#> [[25]]
#> [1] "c" "b" "c"
#> 
#> [[26]]
#> [1] "c" "a" "b"
#> 
#> [[27]]
#> [1] "b" "c" "b"
#> 
#> [[28]]
#> [1] "a" "b" "c"
#> 
#> [[29]]
#> [1] "b" "b" "b"
#> 
#> [[30]]
#> [1] "b" "c" "b"
#> 
#> [[31]]
#> [1] "b" "c" "b"
#> 
#> [[32]]
#> [1] "a" "a" "c"
#> 
#> [[33]]
#> [1] "b" "b" "c"
#> 
#> [[34]]
#> [1] "c" "b" "c"
#> 
#> [[35]]
#> [1] "b" "b" "c"
#> 
#> [[36]]
#> [1] "c" "b" "c"
#> 
#> [[37]]
#> [1] "b" "c" "b"
#> 
#> [[38]]
#> [1] "c" "b" "c"
#> 
#> [[39]]
#> [1] "a" "b" "c"
#> 
#> [[40]]
#> [1] "c" "b" "c"
#> 
#> [[41]]
#> [1] "b" "c" "b"
#> 
#> [[42]]
#> [1] "a" "c" "b"
#> 
#> [[43]]
#> [1] "c" "b" "c"
#> 
#> [[44]]
#> [1] "b" "c" "a"
#> 
#> [[45]]
#> [1] "b" "c" "b"
#> 
#> [[46]]
#> [1] "a" "a" "a"
#> 
#> [[47]]
#> [1] "c" "b" "c"
#> 
#> [[48]]
#> [1] "b" "c" "c"
#> 
#> [[49]]
#> [1] "c" "a" "a"
#> 
#> [[50]]
#> [1] "b" "b" "c"
#> 
#> [[51]]
#> [1] "b" "c" "a"
#> 
#> [[52]]
#> [1] "c" "b" "c"
#> 
#> [[53]]
#> [1] "b" "c" "b"
#> 
#> [[54]]
#> [1] "b" "c" "b"
#> 
#> [[55]]
#> [1] "b" "c" "b"
#> 
#> [[56]]
#> [1] "b" "c" "b"
#> 
#> [[57]]
#> [1] "b" "b" "c"
#> 
#> [[58]]
#> [1] "b" "c" "b"
#> 
#> [[59]]
#> [1] "a" "a" "b"
#> 
#> [[60]]
#> [1] "c" "c" "b"
#> 
#> [[61]]
#> [1] "c" "b" "c"
#> 
#> [[62]]
#> [1] "b" "c" "b"
#> 
#> [[63]]
#> [1] "c" "b" "b"
#> 
#> [[64]]
#> [1] "b" "c" "b"
#> 
#> [[65]]
#> [1] "b" "b" "c"
#> 
#> [[66]]
#> [1] "c" "b" "c"
#> 
#> [[67]]
#> [1] "a" "b" "c"
#> 
#> [[68]]
#> [1] "c" "b" "c"
#> 
#> [[69]]
#> [1] "c" "b" "c"
#> 
#> [[70]]
#> [1] "c" "b" "b"
#> 
#> [[71]]
#> [1] "b" "c" "b"
#> 
#> [[72]]
#> [1] "a" "b" "c"
#> 
#> [[73]]
#> [1] "c" "b" "c"
#> 
#> [[74]]
#> [1] "b" "c" "c"
#> 
#> [[75]]
#> [1] "b" "c" "b"
#> 
#> [[76]]
#> [1] "b" "c" "b"
#> 
#> [[77]]
#> [1] "b" "c" "b"
#> 
#> [[78]]
#> [1] "b" "c" "b"
#> 
#> [[79]]
#> [1] "c" "b" "c"
#> 
#> [[80]]
#> [1] "c" "c" "b"
#> 
#> [[81]]
#> [1] "c" "b" "c"
#> 
#> [[82]]
#> [1] "c" "b" "c"
#> 
#> [[83]]
#> [1] "b" "c" "b"
#> 
#> [[84]]
#> [1] "b" "c" "b"
#> 
#> [[85]]
#> [1] "c" "b" "b"
#> 
#> [[86]]
#> [1] "b" "b" "c"
#> 
#> [[87]]
#> [1] "b" "c" "b"
#> 
#> [[88]]
#> [1] "b" "c" "b"
#> 
#> [[89]]
#> [1] "c" "c" "b"
#> 
#> [[90]]
#> [1] "b" "c" "b"
#> 
#> [[91]]
#> [1] "b" "c" "b"
#> 
#> [[92]]
#> [1] "b" "b" "b"
#> 
#> [[93]]
#> [1] "b" "c" "b"
#> 
#> [[94]]
#> [1] "a" "c" "b"
#> 
#> [[95]]
#> [1] "b" "c" "b"
#> 
#> [[96]]
#> [1] "b" "c" "b"
#> 
#> [[97]]
#> [1] "b" "c" "b"
#> 
#> [[98]]
#> [1] "b" "c" "b"
#> 
#> [[99]]
#> [1] "b" "c" "a"
#> 
#> [[100]]
#> [1] "c" "b" "c"
#> 
     
```
