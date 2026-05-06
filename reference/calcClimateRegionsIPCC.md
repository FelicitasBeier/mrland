# calcClimateRegionsIPCC

calculates IPCC Climate Regions (IPCC2006 ch.4.3) based on t, ppt, pet
from LPJml. elevation dimension not included for tropical montane class

## Usage

``` r
calcClimateRegionsIPCC(
  landusetypes = "all",
  cellular = FALSE,
  yearly = FALSE,
  convert = TRUE
)
```

## Arguments

- landusetypes:

  all or only one

- cellular:

  FALSE for country level, TRUE for cells

- yearly:

  FALSE for normal magpie 5 year time spans, TRUE for yearly

- convert:

  fills missing countries for country level aggregation with warm
  temperate moist (mostly small island nations)

## Value

Country or cellular magpie object with fraction of each climate region
by country or cell

## Author

David Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ClimateRegionsIPCC")
} # }
```
