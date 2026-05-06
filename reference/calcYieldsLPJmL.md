# calcYieldsLPJmL

This function extracts yields from LPJmL for all years

## Usage

``` r
calcYieldsLPJmL(
  lpjml = "ggcmi_phase3_nchecks_bft_e511ac58",
  climatetype = "GSWP3-W5E5:historical",
  cells = "lpjcell"
)
```

## Arguments

- lpjml:

  Defines LPJmL version for main crop inputs

- climatetype:

  Switch between different climate scenarios

- cells:

  if cellular is TRUE: "magpiecell" for 59199 cells or "lpjcell" for
  67420 cells

## Value

magpie object in cellular resolution

## Author

Kristine Karstens, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("YieldsLPJmL", aggregate = FALSE)
} # }
```
