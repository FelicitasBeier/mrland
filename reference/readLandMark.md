# readLandMark

Reads Indigenous Peoples’ and Local Community Lands from LandMark
"LandMark. 2025. Indicative Areas of Indigenous and Community Land
Rights. Data file from LandMark: The Global Platform of Indigenous and
Community Lands. Available at www.landmarkmap.org."

Protected areas (WDPA+China) were masked at a spatial resolution of 10
arcseconds before aggregating the data to 0.5°.

## Usage

``` r
readLandMark(subtype = "delineated")
```

## Arguments

- subtype:

  Whether `"delineated"` or `"indicative"` lands and territories are
  returned.

## Value

Returns magpie objects with the land area covered by the IPLC lands that
was unprotected in 2020.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("LandMark", convert = "onlycorrect")
} # }
```
