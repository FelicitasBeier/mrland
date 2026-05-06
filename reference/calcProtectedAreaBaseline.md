# calcProtectedAreaBaseline

Returns protected land area (Mha) in terms of cropland, pasture, forest
and other land between 1995 and 2020.

## Usage

``` r
calcProtectedAreaBaseline(
  magpie_input = TRUE,
  nclasses = "seven",
  cells = "lpjcell"
)
```

## Arguments

- magpie_input:

  Whether data should be transformed (based on LUH3 data) to match land
  use types used in MAgPIE.

- nclasses:

  If `magpie_input = TRUE`. Options are either "seven" or "nine". Note
  that by default, the protected area is reported for urban land and
  forestry is zero.

  - "seven" separates primary and secondary forest and includes "crop",
    "past", "forestry", "primforest", "secdforest", "urban" and "other"

  - "nine" adds the separation of pasture and rangelands, as well as a
    differentiation of primary and secondary non-forest vegetation and
    therefore returns "crop", "past", "range", "forestry", "primforest",
    "secdforest", "urban", "primother" and "secdother"

- cells:

  (deprecated) only lpjcell (67420 cells)

## Value

List with a magpie object

## See also

[`calcProtectArea`](calcProtectArea.md)

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ProtectedAreaBaseline", aggregate = FALSE)
} # }
```
