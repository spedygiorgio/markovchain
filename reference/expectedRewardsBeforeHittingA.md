# Expected first passage Rewards for a set of states in a markovchain

Given a markovchain object and reward values for every state, function
calculates expected reward value for a set A of states after n steps.

## Usage

``` r
expectedRewardsBeforeHittingA(markovchain, A, state, rewards, n)
```

## Arguments

- markovchain:

  the markovchain-class object

- A:

  set of states for first passage expected reward

- state:

  initial state

- rewards:

  vector depicting rewards coressponding to states

- n:

  no of steps of the process

## Value

returns a expected reward (numerical value) as described above

## Details

The function returns the value of expected first passage rewards given
rewards coressponding to every state, an initial state and number of
steps.

## Author

Sai Bhargav Yalamanchi, Vandit Jain
