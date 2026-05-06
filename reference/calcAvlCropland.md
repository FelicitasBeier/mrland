# calcAvlCropland

Calculates the total available cropland per grid cell, based on physical
cropland suitability data or other criteria, such as constraints on
cropland expansion

## Usage

``` r
calcAvlCropland(
  marginal_land = "magpie",
  cell_upper_bound = 0.9,
  country_level = FALSE,
  cells = "lpjcell",
  luhBaseYear = "y1995"
)
```

## Arguments

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

- cell_upper_bound:

  Upper bound for cropland at the grid cell level. Even if, for
  instance, the total available cropland area equals the land area in a
  grid cell, cropland cannot be expanded above this value.

- country_level:

  Whether output shall be at country level. Requires aggregate=FALSE in
  calcOutput.

- cells:

  (deprecated) only lpjcell (67420 cells)

- luhBaseYear:

  Base year of LUH land area

## Value

magpie object in cellular resolution

## Author

Patrick v. Jeetze, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("AvlCropland", aggregate = FALSE)
} # }
```
