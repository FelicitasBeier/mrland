# calcForestLossByDriver

Forest area lost per year by driver, in Mha, from the Global Forest
Watch driver product of Sims et al. (2025) (`src = "GFW"`) or from Table
1 of Curtis et al. (2018) (`src = "Curtis"`, default).

## Usage

``` r
calcForestLossByDriver(src = "Curtis", period = 2015:2024)
```

## Arguments

- src:

  `"Curtis"` (default) or `"GFW"`

- period:

  Years averaged into the annual rate, GFW only; Curtis carries a single
  2001-2015 mean

## Value

MAgPIE object with forest area lost per year by driver, Mha

## Details

Of the GFW driver classes only shifting cultivation is returned. Curtis
et al. (2018) report loss for seven world regions. It is distributed to
countries by FRA 2020 naturally regenerating forest area in 2010 and
returned with all drivers and with `overall`, the mean annual forest
fire area of FRA 2020: the columns of `f35_forest_lost_share.cs3` that
MAgPIE's `pot_forest_may24` realization declares.

## See also

[`readGFWLossByDriver`](readGFWLossByDriver.md),
[`readForestLossDrivers`](readForestLossDrivers.md)

## Author

Abhijeet Mishra, Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ForestLossByDriver", aggregate = FALSE)
calcOutput("ForestLossByDriver", src = "GFW", aggregate = FALSE)
} # }
```
