# calcMulticroppingYieldIncrease

Calculates yield increase achieved through multiple cropping (as factor
of off season to main season crop yield) under irrigated and rainfed
conditions respectively. Optionally: return which grid cells are
potentially suitable for multiple cropping activities under rainfed and
irrigated conditions. Calculation is based on grassland gross primary
production (GPP) in the growing period of the respective crop and annual
grass GPP.

## Usage

``` r
calcMulticroppingYieldIncrease(
  selectyears,
  lpjml,
  climatetype,
  fallowFactor = 0.75
)
```

## Arguments

- selectyears:

  Years to be returned

- lpjml:

  LPJmL version required for respective inputs as single string: "crop"
  version

- climatetype:

  Switch between different climate scenarios or historical baseline
  "GSWP3-W5E5:historical"

- fallowFactor:

  Factor determining yield reduction in off season due to fallow period
  between harvest of first (main) season and sowing of second (off)
  season

## Value

magpie object in cellular resolution

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("MulticroppingYieldIncrease", aggregate = FALSE)
} # }
```
