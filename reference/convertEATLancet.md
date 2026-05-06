# convertEATLancet

Convert data from the EAT Lancet Commission to be used in MAgPIE

## Usage

``` r
convertEATLancet(x, subtype)
```

## Arguments

- x:

  MAgPIE object containing EAT Lancet data at mixed country-region
  resolution

- subtype:

  Type of EAT Lancet data that should be read. Available types are:

  - `cons_data`: Consumption analysis ("EAT_Lancet_cons_data.csv")

  - `recommend`: Food recommendations ("EAT_Lancet_recommendations.csv")

## Value

EAT Lancet data as MAgPIE object at ISO country level

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource(type = "EATLancet", subtype = "cons_data")
} # }
```
