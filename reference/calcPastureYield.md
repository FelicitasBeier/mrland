# calcPastureYield

Provides pasture yields defined as ratio of grazed biomass to grazed
area

## Usage

``` r
calcPastureYield(range_pastr = FALSE)
```

## Arguments

- range_pastr:

  Boolean value indicating if the grass yields should be split between
  rangelands and pastures.

## Value

Pasture yields and corresponding weights as a list of two MAgPIE objects

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`calcFAOmassbalance`](https://rdrr.io/pkg/mrcommons/man/calcFAOmassbalance.html),
[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl, Marcos Alves

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("PastureYield")
} # }
```
