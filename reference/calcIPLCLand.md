# calcIPLCLand

Returns unprotected land area (Mha) covered by Indigenous Peoples' and
Local Community Lands.

## Usage

``` r
calcIPLCLand(maginput = TRUE, nclasses = "seven", datasource = "LandMark")
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

- datasource:

  Currently only `"LandMark"`

## Value

List with a magpie object

## See also

[`readLandMark`](readLandMark.md)

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("calcIPLCLand", aggregate = FALSE)
} # }
```
