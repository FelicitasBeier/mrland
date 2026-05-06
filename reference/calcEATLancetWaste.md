# calcEATLancetWaste

Calculates the ratio between food supply at household level and food
intake for different MAgPIE commodities accounting for food-specific
estimates of baseline intake of quantification of EAT Lancet diets by
the EAT-Lancet comission, as well as for FAO food waste shares.

## Usage

``` r
calcEATLancetWaste(out_type = "ratio")
```

## Arguments

- out_type:

  ratio: total food supply to total intake. ratio_detailed_calib:
  calibrated food-specific estimates. ratio_detailed: food-specific
  estimates based on FAO food waste shares calib: factor for calibrating
  estimates based on FAO waste shares to food supply

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readEATLancet`](readEATLancet.md),
[`calcEATLancetDiets`](calcEATLancetDiets.md),
[`convertEATLancet`](convertEATLancet.md)

## Author

Isabelle Weindl

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EATLancetWaste")
} # }
```
