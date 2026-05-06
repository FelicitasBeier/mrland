# calcHalfEarth

Function calculates land area in 'Half Earth' conservation priority area

## Usage

``` r
calcHalfEarth(cells = "lpjcell", nclasses = "seven")
```

## Arguments

- cells:

  (deprecated) always lpjcell (67420 cells)

- nclasses:

  Options are either "seven" or "nine".

  - "seven" separates primary and secondary forest and includes "crop",
    "past", "forestry", "primforest", "secdforest", "urban" and "other"

  - "nine" adds the separation of pasture and rangelands, as well as a
    differentiation of primary and secondary non-forest vegetation and
    therefore returns "crop", "past", "range", "forestry", "primforest",
    "secdforest", "urban", "primother" and "secdother"

## Value

magpie object in cellular resolution with different protection options
in conservation priority areas

## Author

Patrick v. Jeetze, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("HalfEarth", aggregate = FALSE)
} # }
```
