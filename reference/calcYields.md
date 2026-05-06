# calcYields

This function extracts yields from LPJmL and transforms them to MAgPIE
crops calibrating proxy crops to FAO yields. Optionally, ISIMIP yields
can be returned.

## Usage

``` r
calcYields(
  source = c(lpjml = "ggcmi_phase3_nchecks_9ca735cb", isimip = NULL),
  climatetype = "GSWP3-W5E5:historical",
  cells = "lpjcell",
  selectyears = seq(1965, 2100, by = 5),
  weighting = "totalCrop",
  multicropping = FALSE,
  indiaYields = FALSE,
  scaleFactor = 0.3,
  marginal_land = "magpie"
)
```

## Arguments

- source:

  Defines LPJmL version for main crop inputs and isimip replacement. For
  isimip choose crop model/gcm/rcp/co2 combination formatted like this:
  "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b"

- climatetype:

  Switch between different climate scenarios

- cells:

  if cellular is TRUE: "magpiecell" for 59199 cells or "lpjcell" for
  67420 cells

- selectyears:

  Years to be returned

- weighting:

  use of different weights (totalCrop (default), totalLUspecific,
  cropSpecific, crop+irrigSpecific, avlCropland, avlCropland+avlPasture)

- multicropping:

  Multicropping activated (TRUE) or not (FALSE) and Multiple Cropping
  Suitability mask selected (mask can be: "none": no mask applied (only
  for development purposes) "actual:total": currently multicropped areas
  calculated from total harvested areas and total physical areas per
  cell from readLandInG "actual:crop" (crop-specific),
  "actual:irrigation" (irrigation-specific), "actual:irrig_crop" (crop-
  and irrigation-specific), "potential:endogenous": potentially
  multicropped areas given temperature and productivity limits
  "potential:exogenous": potentially multicropped areas given GAEZ
  suitability classification) (e.g. TRUE:actual:total; TRUE:none; FALSE)

- indiaYields:

  if TRUE returns scaled yields for rainfed crops in India

- scaleFactor:

  integer value by which indiaYields will be scaled

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
calcOutput("Yields", aggregate = FALSE)
} # }
```
