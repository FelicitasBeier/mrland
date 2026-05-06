# calcPhotosynthesisTemperature

This function calculates crop-specific temperature limits for the
multicropping mask based on the photosynthesis optimum and the LPJmL
parameters temp_co2 and temp_photos

## Usage

``` r
calcPhotosynthesisTemperature(threshold = 0.8)
```

## Arguments

- threshold:

  Photosynthesis efficiency threshold (between 0 and 1)

## Value

magpie object

## Author

Felicitas Beier, Jens Heinke

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("PhotosynthesisTemperature", aggregate = FALSE)
} # }
```
