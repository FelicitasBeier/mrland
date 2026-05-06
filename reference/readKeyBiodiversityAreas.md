# readKeyBiodiversityAreas

Reads land area covered by for Key Biodiversity Areas
(https://www.keybiodiversityareas.org/) that was unprotected in 2020.
Protected areas were masked at a spatial resolution of 10 arc seconds
before aggregating the data to 0.5°.

## Usage

``` r
readKeyBiodiversityAreas(subtype = "unprotected")
```

## Arguments

- subtype:

  "unprotected" or "all"

## Value

Returns magpie objects with the area covered by unprotected Key
Biodiversity Areas per grid cell

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("KeyBiodiversityAreas", convert = "onlycorrect")
} # }
```
