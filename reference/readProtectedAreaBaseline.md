# readProtectedAreaBaseline

Reads spatial land cover information within protected areas. Land cover
information for protected areas has been extracted from ESA CCI land
use/land cover data (https://www.esa-landcover-cci.org/) and data from
the WDPA data base (https://www.protectedplanet.net).

## Usage

``` r
readProtectedAreaBaseline()
```

## Value

Returns magpie object with the protected area separated for each land
type (cropland, pasture, forest, other land) per grid cell from 1995 to
2020.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("ProtectedAreaBaseline", convert = "onlycorrect")
} # }
```
