# readNoon2022

Reads irrecoverble carbon data set published by Noon, M. L., Goldstein,
A., Ledezma, J. C., Roehrdanz, P. R., Cook-Patton, S. C., Spawn-Lee, S.
A., Wright, T. M., Gonzalez-Roglich, M., Hole, D. G., Rockström, J., &
Turner, W. R. (2022). Mapping the irrecoverable carbon in Earth’s
ecosystems. Nature Sustainability, 5(1), Article 1.
https://doi.org/10.1038/s41893-021-00803-6 Protected areas were masked
at a spatial resolution of 10 arc seconds before aggregating the data to
0.5°.

## Usage

``` r
readNoon2022(subtype = "land:IrrC_50pc")
```

## Arguments

- subtype:

  Defines whether carbon data or land area and related subtypes should
  be returned (see options below). Carbon or land subtypes need to be
  specified via ":" The different subtypes for land are: "IrrC_30pc",
  "IrrC_40pc", "IrrC_50pc", "IrrC_60pc", "IrrC_70pc", "IrrC_80pc",
  "IrrC_90pc", "IrrC_100pc" which corresponds to the land area that was
  unprotected in 2020 and is covered by the respective percentile of all
  irrecoverable carbon. `IrrC_50pc` e.g. returns all unprotected land
  that contains the top 50\\% of global irrecoverable carbon.

## Value

Returns magpie objects with the area of unprotected irrecoverable carbon
land per grid cell

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Noon2022", convert = "onlycorrect")
} # }
```
