# readForestLossDrivers

Reads Table 1 of Curtis et al. (2018), tree cover loss 2001-2015 by
world region and dominant driver, transcribed by hand from the article
into `forest_loss_fixed.csv`.

## Usage

``` r
readForestLossDrivers()
```

## Value

magpie object with the mean annual loss per Curtis region and driver,
Mha

## References

Curtis, P. G., Slay, C. M., Harris, N. L., Tyukavina, A. and Hansen, M.
C. (2018) Classifying drivers of global forest loss. Science 361,
1108-1111. doi:10.1126/science.aau3445

## See also

[`calcForestLossByDriver`](calcForestLossByDriver.md)

## Author

Abhijeet Mishra, Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("ForestLossDrivers", convert = FALSE)
} # }
```
