# calcBHIFL

Function calculates land area in conservation priority areas

## Usage

``` r
calcBHIFL(cells = "lpjcell", nclasses = "seven")
```

## Arguments

- cells:

  number of cells of landmask (select "magpiecell" for 59199 cells or
  "lpjcell" for 67420 cells)

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
calcOutput("BHIFL", aggregate = FALSE)
} # }
```
