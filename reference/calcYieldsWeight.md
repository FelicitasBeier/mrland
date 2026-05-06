# calcYieldsWeight

This function calculates the crop area weightings to use for yields.

## Usage

``` r
calcYieldsWeight(
  cells = "lpjcell",
  weighting = "totalCrop",
  marginal_land = "magpie"
)
```

## Arguments

- cells:

  if cellular is TRUE: "magpiecell" for 59199 cells or "lpjcell" for
  67420 cells

- weighting:

  use of different weights (totalCrop (default), totalLUspecific,
  cropSpecific, crop+irrigSpecific, avlCropland, avlCropland+avlPasture)

- marginal_land:

  Defines which share of marginal land should be included (see options
  below) and whether suitable land under irrigated conditions
  ("irrigated"), under rainfed conditions ("rainfed") or suitability
  under rainfed conditions including currently irrigated land
  (rainfed_and_irrigated) should be used. Options combined via ":" The
  different marginal land options are:

  - `"all_marginal"`: All marginal land (suitability index between
    0-0.33) is included as suitable

  - `"q33_marginal"`: The bottom tertile (suitability index below 0.13)
    of the marginal land area is excluded.

  - `"q50_marginal"`: The bottom half (suitability index below 0.18) of
    the marginal land area is excluded.

  - `"q66_marginal"`: The first and second tertile (suitability index
    below 0.23) of the marginal land area are excluded.

  - `"q75_marginal"`: The first, second and third quartiles (suitability
    index below 0.25) of the marginal land are are excluded

  - `"no_marginal"`: Areas with a suitability index of 0.33 and lower
    are excluded.

  - `"magpie"`: Returns "all_marginal:rainfed_and_irrigated",
    "q33_marginal:rainfed_and_irrigated" and
    "no_marginal:rainfed_and_irrigated" in a magclass object to be used
    as magpie input.

## Value

magpie object in cellular resolution

## Author

Kristine Karstens, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("YieldsWeight", yields, aggregate = FALSE)
} # }
```
