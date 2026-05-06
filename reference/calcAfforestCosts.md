# Aggregation and calculation of the mean of each MAgPIE region for the source SathayeForest

This function aggregates the data from source SathayeForest. A weight is
implemented as the mean for each MAgPIE region is calculated.

## Usage

``` r
calcAfforestCosts()
```

## Value

MAgPIE object of the calculated means of each MAgPIE region

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readSathayeForest`](readSathayeForest.md),
[`convertSathayeForest`](convertSathayeForest.md)

## Author

Nele Steinmetz

## Examples

``` r

if (FALSE) { # \dontrun{

a <- calcOutput("AfforestCosts")

} # }
```
