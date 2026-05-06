# calcEATFruitvegRatio

Calculates the share of fruits and vegetables in the calorie supply from
the others MAgPIE commodity for the past. Information on the calorie
supply from fruits and vegetables is relevant in the context of dietary
recommendations, e.g. as proposed by the EAT-Lancet Commission on
healthy diets from sustainable food systems.

## Usage

``` r
calcEATFruitvegRatio(populationweight = "PopulationPast")
```

## Arguments

- populationweight:

  datasource of populationweight: FAO can be selected in order to better
  meet exact values. Normal datasource is PopulationPast

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`calcEATLancetTargets`](calcEATLancetTargets.md),
[`calcFAOharmonized`](https://rdrr.io/pkg/mrfaocore/man/calcFAOharmonized.html),
[`calcEATLancetDiets`](calcEATLancetDiets.md)

## Author

Isabelle Weindl, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EATFruitvegRatio")
} # }
```
