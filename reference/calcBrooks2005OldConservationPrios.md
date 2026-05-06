# calcBrooks2005OldConservationPrios

Function calculates land area in conservation priority areas

## Usage

``` r
calcBrooks2005OldConservationPrios(cells = "lpjcell", nclasses = "seven")
```

## Arguments

- cells:

  (deprecated) number of cells of landmask ("lpjcell" for 67420 cells)

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

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ConservationPriority", aggregate = FALSE)
} # }
```
