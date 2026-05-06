# calcAtmosphericDepositionRates

Conputes Atmospheric (nitrogen) deposition rates per area on different
land-use types.

## Usage

``` r
calcAtmosphericDepositionRates(cellular = FALSE, cells = "lpjcell")
```

## Arguments

- cellular:

  TRUE for results on 0.5 degree grid.

- cells:

  magpiecell (59199 cells) or lpjcell (67420 cells)

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcAtmosphericDeposition`](https://rdrr.io/pkg/mrcommons/man/calcAtmosphericDeposition.html),
[`calcNitrogenBudgetCropland`](https://rdrr.io/pkg/mrcommons/man/calcNitrogenBudgetCropland.html)

## Author

Benjamin Leon Bodirsky

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AtmosphericDepositionRates")
} # }
```
