# convertAQUASTAT

Convert data based on AQUASTAT database
(http://www.fao.org/nr/water/aquastat/data/query/index.html?lang=en)

## Usage

``` r
convertAQUASTAT(x, subtype)
```

## Arguments

- x:

  MAgPIE object containing AQUASTAT data on country level

- subtype:

  - `ConsAgri`: 4454\|Conservation agriculture area (1000 ha)
    4454_conservation_agriculture_area_in_1000_ha.csv

  - `ConsAgriShare`: 4455\|Commoditiy Balance LivestockConservation
    agriculture area as 4455_conservation_agriculture_area_as_share_of\_
    arable_land_areas.csv)

  - `rf2irRatio`: Ratio between rainfed and irrigated yields (
    Ratio_between_rainfed_and_irrigated_yields.csv

## Value

magpie objects with results on contury level

## Author

Kristine Karstens

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("AQUASTAT", subtype = "ConsAgri", convert = TRUE)
} # }
```
