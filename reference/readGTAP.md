# readGTAP

Read BaseData and BaseView in GTAP database that has been downlodaded
from the GTAP wewbsite.

## Usage

``` r
readGTAP(subtype = NULL)
```

## Arguments

- subtype:

  Type of GTAP data that should be read. So far available are:

  - GTAP7:

    - `GTAP7_VIWS`: Trade - Bilateral Imports at World Prices

    - `GTAP7_VIMS`: Trade - Bilateral Imports at Market Prices

    - `GTAP7_VXWD`: Trade - Bilateral Exports at World Prices

    - `GTAP7_VXMD`: Trade - Bilateral Exports at Market Prices

    - `GTAP7_VDFM`: Intermediates - Firms' Domestic Purchases at Market
      Prices

    - `GTAP7_VIFM`: Intermediates - Firms' Imports at Market Prices

    - `GTAP7_VFM`: Endowments - Firms' Purchases at Market Prices

    - `GTAP7_VOA`: Payment received by producers (fram gtate value)

    - `GTAP7_VOM`: Value of output at dometic market prices

  - GTAP8:

    - `GTAP8_VIWS`: Trade - Bilateral Imports at World Prices

    - `GTAP8_VIMS`: Trade - Bilateral Imports at Market Prices

    - `GTAP8_VXWD`: Trade - Bilateral Exports at World Prices

    - `GTAP8_VXMD`: Trade - Bilateral Exports at Market Prices

    - `GTAP8_VDFM`: Intermediates - Firms' Domestic Purchases at Market
      Prices

    - `GTAP8_VIFM`: Intermediates - Firms' Imports at Market Prices

    - `GTAP8_VFM`: Endowments - Firms' Purchases at Market Prices

    - `GTAP8_VOA`: Payment received by producers (fram gtate value)

    - `GTAP8_VOM`: Value of output at dometic market prices

## Value

GTAP data as a MAgPie-Object

## Author

Stephen Wirth, Xiaoxi Wang

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("GTAP7", "VIWS")
} # }
```
