# calcISIMIP3bYields

reads and cleans up ISIMIP3b crop yield data

## Usage

``` r
calcISIMIP3bYields(
  subtype = "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b",
  smooth = TRUE,
  cells = "lpjcell"
)
```

## Arguments

- subtype:

  subtype of yield based on readISIMIPoutputs, for crop yields

- smooth:

  smooth cells via spline

- cells:

  magpie or lpjcell

## Value

magpie object in cellular resolution

## Author

David Meng-Chuen Chen, Edna Molina Bacca

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ISIMIP3bYields", aggregate = FALSE)
} # }
```
