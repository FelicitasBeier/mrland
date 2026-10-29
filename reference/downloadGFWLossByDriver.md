# downloadGFWLossByDriver

Downloads tree cover loss by country, year and dominant driver from the
Global Forest Watch data API, table `gadm__tcl__iso_change`: the
UMD/Hansen annual tree cover loss cross-tabulated with the 1 km
dominant-driver map of Sims et al. (2025).

## Usage

``` r
downloadGFWLossByDriver()
```

## See also

[`readGFWLossByDriver`](readGFWLossByDriver.md)

## Author

Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
madrat::downloadSource("GFWLossByDriver")
} # }
```
