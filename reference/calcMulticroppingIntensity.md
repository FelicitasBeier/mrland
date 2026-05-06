# calcMulticroppingIntensity

Returns cropping intensity according to LandInG data given the chosen
scenario

## Usage

``` r
calcMulticroppingIntensity(scenario, selectyears, sectoral = "lpj")
```

## Arguments

- scenario:

  "total": currently multicropped areas calculated from total harvested
  areas and total physical areas per cell from readLandInG "crop"
  (crop-specific), "irrigation" (irrigation-specific), "irrig_crop"
  (crop- and irrigation-specific)

- selectyears:

  Years to be returned

- sectoral:

  "kcr" MAgPIE crop types, and "lpj" LPJmL crop types

## Value

magpie object in cellular resolution

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("MulticroppingIntensity", aggregate = FALSE)
} # }
```
