# Convert data from the NIN Lancet Comission

Convert data from the NIN Lancet Comission to ISO country level.

## Usage

``` r
convertNIN(x, subtype)
```

## Arguments

- x:

  MAgPIE object containing NIN Lancet data at mixed country-region
  resolution

- subtype:

  Type of NIN Lancet data that should be read. Available types are:

  - `cons_data`: Consumption analysis ("NIN_Lancet_cons_data.csv")

  - `recommend`: Food recommendations ("NIN_recommendations.csv")

## Value

NIN Lancet data as MAgPIE object at ISO country level

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl

## Examples

``` r

if (FALSE) { # \dontrun{
a <- readSource(type="NIN",subtype="cons_data")
} # }
```
