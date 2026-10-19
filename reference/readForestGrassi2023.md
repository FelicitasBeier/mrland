# readForestGrassi2023

Reads the 0.5 degree map of intact and non-intact forest area of Grassi
et al. (2023), Harmonising the land-use flux estimates of global models
and national inventories for 2000-2020, Earth Syst. Sci. Data 15,
1093-1114, https://doi.org/10.5194/essd-15-1093-2023. Intact forest
follows Potapov et al. (2017) for 2013, with the masks for Canada,
Brazil and Russia updated from national information; non-intact forest
is forest after Hansen et al. (2013) minus intact forest.

## Usage

``` r
readForestGrassi2023()
```

## Value

magpie object with intact and non-intact forest area (Mha) on the 67420
cells

## See also

[`downloadForestGrassi2023`](downloadForestGrassi2023.md),
[`correctForestGrassi2023`](correctForestGrassi2023.md)

## Author

Florian Humpenöder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("ForestGrassi2023", convert = "onlycorrect")
} # }
```
