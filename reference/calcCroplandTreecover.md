# calcCroplandTreecover

Returns area on cropland covered by trees (Mha).

## Usage

``` r
calcCroplandTreecover(maginput = TRUE, countryLevel = FALSE)
```

## Arguments

- maginput:

  Whether data should be corrected to align with cropland initialised in
  MAgPIE.

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
