# calcEATLancetTargets

Calculates minimum and maximum targets for healthy food intake according
to reference recommendations proposed by the EAT-Lancet Commission on
healthy diets from sustainable food systems, specified for different
MAgPIE commodities.

## Usage

``` r
calcEATLancetTargets(attributes = "kcal/d")
```

## Arguments

- attributes:

  Attributes of food commodities (available: kcal/d and g/d)

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readEATLancet`](readEATLancet.md),
[`calcEATLancetDiets`](calcEATLancetDiets.md)

## Author

Isabelle Weindl

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EATLancetTargets")
} # }
```
