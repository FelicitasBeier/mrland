# calcNINDiets

Calculates daily per capita intake for MAgPIE food commodities that are
consistent with diet scenarios developed by the NIN-Lancet Commission on
healthy diets from sustainable food systems. The unit is kcal/day per
capita or wm/day per capita. Mapping of intake from NIN Lancet to MAgPIE
food commodities is done indivudually for the different available units.

## Usage

``` r
calcNINDiets(attributes = c("wm", "kcal"), calib = TRUE, FAOcountr = FALSE)
```

## Arguments

- attributes:

  attributes of different food commodities (available: kcal and wm).

- calib:

  if TRUE, total daily per capita intake for MAgPIE food commodities is
  calibrated to NIN Lancet total intake.

- FAOcountr:

  if TRUE, estimates for countries not covered in FAOSTAT are set to
  Zero.

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readNIN`](readNIN.md), [`convertNIN`](convertNIN.md)

## Author

Isabelle Weindl

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("NINDiets")
} # }
```
