# calcCroplandTreecover

Returns area on cropland covered by trees (Mha).

## Usage

``` r
calcCroplandTreecover(
  maginput = TRUE,
  cells = "magpiecell",
  countryLevel = FALSE
)
```

## Arguments

- maginput:

  Whether data should be corrected to align with cropland initialised in
  MAgPIE.

- cells:

  magpiecell (59199 cells) or lpjcell (67420 cells)

- countryLevel:

  Whether output shall be at country level. Requires aggregate=FALSE in
  calcOutput.

## Value

List with a magpie object

## See also

[`readCopernicus`](readCopernicus.md)

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("CroplandTreecover", aggregate = FALSE)
} # }
```
