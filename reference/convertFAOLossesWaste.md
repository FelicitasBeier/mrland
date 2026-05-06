# Convert data on food losses and waste from FAO for several commodity groups

Convert data on food losses and waste on ISO country level.

## Usage

``` r
convertFAOLossesWaste(x, subtype)
```

## Arguments

- x:

  MAgPIE object containing data on food losses and waste at mixed
  country-region resolution

- subtype:

  Steps of the food supply chain where food losses and waste occur.
  Available types are:

  - `Consumption`: consumption level

## Value

Data on food losses and waste as MAgPIE object at ISO country level

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl

## Examples

``` r

if (FALSE) { # \dontrun{
a <- readSource(type="FAOLossesWaste",subtype="Consumption")
} # }
```
