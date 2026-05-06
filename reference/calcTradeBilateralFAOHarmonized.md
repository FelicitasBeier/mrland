# calcTradeBilateralFAOHarmonized

Harmonizes bilateral trade values to FAO mass balance imports often
termed "Import Dependency Ratio"

## Usage

``` r
calcTradeBilateralFAOHarmonized(yearly = TRUE)
```

## Arguments

- yearly:

  whether to calculate yearly data or only magpie 5year timesteps

## Value

Self import to dupply ratio

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`calcFAOmassbalance`](https://rdrr.io/pkg/mrcommons/man/calcFAOmassbalance.html)

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
a <- calcTradeBilateralFAOHarmonized()
} # }
```
