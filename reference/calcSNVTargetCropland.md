# calcSNVTargetCropland

Returns cropland area (Mha) that requires relocation in response of
maintaining 20

## Usage

``` r
calcSNVTargetCropland(maginput = TRUE, cells = "magpiecell")
```

## Arguments

- maginput:

  Whether data should be corrected to align with cropland initialised in
  MAgPIE.

- cells:

  magpiecell (59199 cells) or lpjcell (67420 cells)

## Value

List with a magpie object

## See also

[`readCopernicus`](readCopernicus.md)

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("SNVTargetCropland", aggregate = FALSE)
} # }
```
