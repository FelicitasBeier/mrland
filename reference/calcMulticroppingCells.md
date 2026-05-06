# calcMulticroppingCells

Returns grid cells and crops where multiple cropping takes place given
the chosen scenario

## Usage

``` r
calcMulticroppingCells(
  selectyears,
  lpjml,
  climatetype,
  scenario,
  sectoral = "kcr"
)
```

## Arguments

- selectyears:

  Years to be returned

- lpjml:

  LPJmL version required for respective inputs: natveg or crop

- climatetype:

  Switch between different climate scenarios or historical baseline
  "GSWP3-W5E5:historical"

- scenario:

  "actual:total": currently multicropped areas calculated from total
  harvested areas and total physical areas per cell from readLandInG
  "actual:crop" (crop-specific), "actual:irrigation"
  (irrigation-specific), "actual:irrig_crop" (crop- and
  irrigation-specific) "total" "potential:endogenous": potentially
  multicropped areas given temperature and productivity limits
  "potential:exogenous": potentially multicropped areas given GAEZ
  suitability classification

- sectoral:

  "kcr" MAgPIE crops, and "lpj" LPJmL crops

## Value

magpie object in cellular resolution

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("MulticroppingCells", aggregate = FALSE)
} # }
```
