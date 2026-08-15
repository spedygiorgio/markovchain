# Calculating probability from a ctmc object

This function returns the probability of every state at time t under
different conditions

## Usage

``` r
probabilityatT(C,t,x0,useRCpp)
```

## Arguments

- C:

  A CTMC S4 object

- t:

  final time t

- x0:

  initial state

- useRCpp:

  logical whether to use RCpp implementation

## Value

returns a vector or a matrix in case `x0` is provided or not
respectively.

## Details

The initial state is not mandatory, In case it is not provided, function
returns a matrix of transition function at time `t` else it returns
vector of probaabilities of transition to different states if initial
state was `x0`

## References

INTRODUCTION TO STOCHASTIC PROCESSES WITH R, ROBERT P. DOBROW, Wiley

## Author

Vandit Jain

## Examples

``` r
states <- c("a","b","c","d")
byRow <- TRUE
gen <- matrix(data = c(-1, 1/2, 1/2, 0, 1/4, -1/2, 0, 1/4, 1/6, 0, -1/3, 1/6, 0, 0, 0, 0),
nrow = 4,byrow = byRow, dimnames = list(states,states))
ctmc <- new("ctmc",states = states, byrow = byRow, generator = gen, name = "testctmc")
probabilityatT(ctmc,1,useRCpp = TRUE)
#>            a          b         c          d
#> a 0.41546882 0.24714119 0.2703605 0.06702946
#> b 0.12357060 0.63939068 0.0348290 0.20220972
#> c 0.09012017 0.02321933 0.7411205 0.14553997
#> d 0.00000000 0.00000000 0.0000000 1.00000000
```
