# calcFAOLossesWaste

Calculates the ratio between food supply at household level and food
intake for different MAgPIE commodities based on estimated/assumed FAO
waste shares for each commodity group (optionally also including food
conversion factors into edible matter).

## Usage

``` r
calcFAOLossesWaste(out_type = "waste")
```

## Arguments

- out_type:

  waste: food-specific ratios based on FAO food waste shares
  waste_edible: food-specific ratios based on FAO food waste shares
  including conversion into edible matter

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readFAOLossesWaste`](readFAOLossesWaste.md),
[`calcEATLancetWaste`](calcEATLancetWaste.md)

## Author

Isabelle Weindl

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("FAOLossesWaste")
} # }
```
