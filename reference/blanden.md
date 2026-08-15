# Mobility between income quartiles

This table show mobility between income quartiles for father and sons
for the 1970 cohort born

## Usage

``` r
data(blanden)
```

## Format

An object of class `table` with 4 rows and 4 columns.

## Source

Personal reworking

## Details

The rows represent fathers' income quartile when the son is aged 16,
whilst the columns represent sons' income quartiles when he is aged 30
(in 2000).

## References

Jo Blanden, Paul Gregg and Stephen Machin, Intergenerational Mobility in
Europe and North America, Center for Economic Performances (2005)

## Examples

``` r
data(blanden)
mobilityMc<-as(blanden, "markovchain")
```
