# correctForestGrassi2023

Sets cells without data in the forest map of Grassi et al. (2023) to
zero.

## Usage

``` r
correctForestGrassi2023(x)
```

## Arguments

- x:

  magpie object provided by the read function

## Value

magpie object on cellular level

## See also

[`readForestGrassi2023`](readForestGrassi2023.md)

## Author

Florian Humpenöder

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("ForestGrassi2023", convert = "onlycorrect")
} # }
```
