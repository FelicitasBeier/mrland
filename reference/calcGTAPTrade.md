# calcGTAPTrade

calculate trade data from GTAP dataset

## Usage

``` r
calcGTAPTrade(subtype = NULL, bilateral = FALSE)
```

## Arguments

- subtype:

  GTAP version and subtype, separated by "\_" available versions are
  "GTAP7", "GTAP8", and "GTAP9" GTAP sheets relevant for trade are
  "VIWS": Trade - Bilateral Imports at World Prices "VIMS": Trade -
  Bilateral Imports at Market Prices "VXWD": Trade - Bilateral Exports
  at World Prices "VXMD": Trade - Bilateral Exports at Market Prices
  "VDFM": Intermediates - Firms' Domestic Purchases at Market Prices
  "VIFM": Intermediates - Firms' Imports at Market Prices "VFM":
  Endowments - Firms' Purchases at Market Prices "VOA": Payment received
  by producers (fram gtate value) "VOM": Value of output at dometic
  market prices

- bilateral:

  whether bilateral trade data should be calculated

## Value

Trade related data as an MAgPIE object

## Author

Xiaoxi Wang, David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
x <- calcGTAP("GTAP7_VXMD")
} # }
```
