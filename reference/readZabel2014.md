# readZabel2014

Reads crop suitability data published in Zabel, F., Putzenlechner, B., &
Mauser, W. (2014). Global Agricultural Land Resources – A High
Resolution Suitability Evaluation and Its Perspectives until 2100 under
Climate Change Conditions. PLOS ONE, 9(9), e107522.
https://doi.org/10.1371/journal.pone.0107522 and extracts the share of
suitable cropland per grid cell, depending on different suitability
thresholds.

## Usage

``` r
readZabel2014(subtype = "all_marginal:rainfed_and_irrigated")
```

## Arguments

- subtype:

  Defines which share of marginal land should be included (see options
  below) and whether suitable land under irrigated conditions
  ("irrigated"), under rainfed conditions ("rainfed") or suitability
  under rainfed conditions including currently irrigated land
  (rainfed_and_irrigated) should be used. Options combined via ":" The
  different marginal land options are:

  - `"all_marginal"`: All marginal land (suitability index between
    0-0.33) is included as suitable

  - `"q33_marginal"`: The bottom tertile (suitability index below 0.13)
    of the marginal land () area is excluded.

  - `"q50_marginal"`: The bottom half (suitability index below 0.18) of
    the marginal land area is excluded.

  - `"q66_marginal"`: The first and second tertile (suitability index
    below 0.23) of the marginal land area are excluded.

  - `"q75_marginal"`: The first, second and third quartiles (suitability
    index below 0.25) of the marginal land are are excluded

  - `"no_marginal"`: Areas with a suitability index of 0.33 and lower
    are excluded.

## Value

Returns magpie objects with the share of suitable cropland per grid cell

## Author

Patrick v. Jeetze, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Zabel2014", subtype = "all_marginal:rainfed_and_irrigated", convert = "onlycorrect")
} # }
```
