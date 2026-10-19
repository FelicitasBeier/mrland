# calcForestGrassi2023

Intact and non-intact forest area at 0.5 degree after Grassi et al.
(2023). Intact forest follows Potapov et al. (2017) for 2013, with the
masks for Canada, Brazil and Russia updated from national information.
Grassi et al. use non-intact forest as a proxy of the managed forest
that national greenhouse gas inventories report.

## Usage

``` r
calcForestGrassi2023()
```

## Value

List with a magpie object of intact and non-intact forest area (Mha) on
cellular level

## See also

[`readForestGrassi2023`](readForestGrassi2023.md)

## Author

Florian Humpenöder

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ForestGrassi2023", aggregate = FALSE)
} # }
```
