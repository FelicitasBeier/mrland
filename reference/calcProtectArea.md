# calcProtectArea

Function extracts conservation protected area

## Usage

``` r
calcProtectArea(cells = "lpjcell", bhifl = TRUE)
```

## Arguments

- cells:

  number of cells of landmask (select "magpiecell" for 59199 cells or
  "lpjcell" for 67420 cells)

- bhifl:

  should be TRUE (including BH_IFL scenario) for cellular preprocessing
  revisions \> 4.65

## Value

magpie object in cellular resolution with different protection scenarios

## Author

Felicitas Beier, David Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ProtectArea", aggregate = FALSE)
} # }
```
