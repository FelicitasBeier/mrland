# calcProcessing_balanceflow

Calculates the difference between what is recorded in FAOSTAT production
of secondary products, and what is calculated from our global processing
conversion factors in order to meet FAOSTAT production. Fades out by
2050.

## Usage

``` r
calcProcessing_balanceflow()
```

## Value

magpie object of secondary product balanceflow

## See also

[`calcProcessing_conversion_factors`](calcProcessing_conversion_factors.md)

## Author

Benjamin Bodirsky, David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
a <- calcOutput("Processing_balanceflow",
                aggregate = FALSE)
} # }
```
