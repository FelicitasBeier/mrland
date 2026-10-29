# calcForestLossShare

Calculates which share of forest area is lost per year to each driver of
forest loss.

## Usage

``` r
calcForestLossShare(src = "Curtis", period = 2015:2024)
```

## Arguments

- src:

  Driver data to use, `"Curtis"` (default) or `"GFW"`, see
  [`calcForestLossByDriver`](calcForestLossByDriver.md)

- period:

  Years averaged into the annual rate, GFW only

## Value

MAgPIE object with the share of forest area lost per year by driver

## Details

The denominator is FRA 2020 naturally regenerating forest area, the
natural forest MAgPIE applies the share to: 2020 for `src = "GFW"`, 2010
for `src = "Curtis"` (loss 2001-2015). Shares above 1 are clamped to 1.

## See also

[`calcForestLossByDriver`](calcForestLossByDriver.md)

## Author

Abhijeet Mishra, Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ForestLossShare", aggregate = FALSE)
calcOutput("ForestLossShare", src = "GFW", aggregate = FALSE)
} # }
```
