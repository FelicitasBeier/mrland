# readEATLancet

Read in data from the EAT-Lancet Commission

Read in data from: Food in the Anthropocene: the EAT-Lancet Commission
on healthy diets from sustainable food systems, Lancet 2019
https://doi.org/10.1016/S0140-6736(18)31788-4

## Usage

``` r
readEATLancet(subtype)
```

## Arguments

- subtype:

  Type of EAT-Lancet data that should be read. Available types are:

  - `cons_data`: Consumption analysis ("EAT_Lancet_cons_data.csv")

  - `recommend`: Food recommendations ("EAT_Lancet_recommendations.csv")

## Value

magpie object containing EAT-Lancet Comission data

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl, Jan Philipp Dietrich, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource(type = "EATLancet", subtype = "cons_data")
} # }
```
