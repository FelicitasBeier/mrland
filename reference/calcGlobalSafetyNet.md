# calcGlobalSafetyNet

Returns unprotected land area (Mha) within the Global Safety Net
(Dinerstein et al. 2020).

## Usage

``` r
calcGlobalSafetyNet(maginput = TRUE, nclasses = "seven", cells = "lpjcell")
```

## Arguments

- maginput:

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

  (deprecated) always lpjcell (67420 cells)

## Value

List with a magpie object

## See also

[`readDinerstein2020`](readDinerstein2020.md)

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("calcGlobalSafetyNet", aggregate = FALSE)
} # }
```
