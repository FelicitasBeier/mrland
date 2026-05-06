# readDinerstein2020

Reads Global Safety Net data set published by Dinerstein, E., Joshi, A.
R., Vynne, C., Lee, A. T. L., Pharand-Deschênes, F., França, M.,
Fernando, S., Birch, T., Burkart, K., Asner, G. P., & Olson, D. (2020).
A “Global Safety Net” to reverse biodiversity loss and stabilize Earth’s
climate. Science Advances, 6(36), eabb2824.
https://doi.org/10.1126/sciadv.abb2824

Protected areas and Key Biodiversity Areas were masked at a spatial
resolution of 10 arc seconds before aggregating the data to 0.5°.

## Usage

``` r
readDinerstein2020(subtype = "GSN:distinct_species_assemblages")
```

## Arguments

- subtype:

  Defines which cluster (see Dinerstein et al. 2020) of the Global
  Safety Net is returned.The different subtypes for land are:
  `"GSN:distinct_species_assemblages"`, `"GSN:rare_phenomena"`,
  `"GSN:areas_of_intactness"`, `"GSN:climate_stabilisation_tier1"` and
  `"GSN:climate_stabilisation_tier2"`.

## Value

Returns magpie objects with the land area covered by the Global Safety
Net that is NOT already covered by Key Biodiversity Areas and was
unprotected in 2020.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Dinerstein2020", convert = "onlycorrect")
} # }
```
