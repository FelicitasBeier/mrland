# readCopernicus

Reads either information on the area on cropland covered by trees or
information the cropland area that requires relocation in response of
increasing semi-natural vegetation in farmed landscapes. The data was
derived from high resolution land cover information (LC100) from the
Copernicus Global Land Service. (https://zenodo.org/records/3939050)

## Usage

``` r
readCopernicus(subtype = "CroplandTreecover")
```

## Arguments

- subtype:

  For cropland area covered by trees choose `"CroplandTreecover"`. For
  cropland area requiring relocation in response to increasing SNV
  choose `"SNVTargetCropland"`.

## Value

Returns magpie objects with cropland area covered by trees or cropland
area requiring relocation in order to increase SNV in farmed landscapes.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Copernicus", subtype = "CroplandTreecover", convert = "onlycorrect")
} # }
```
