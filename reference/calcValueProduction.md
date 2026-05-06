# calcValueProduction

calculates production value based on production and prices, only works
for FAO dataset currently

## Usage

``` r
calcValueProduction(datasource = "FAO", cellular = TRUE)
```

## Arguments

- datasource:

  Options of the source of the price data: only FAO has country level
  data

- cellular:

  cellular or iso country values

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcProduction`](https://rdrr.io/pkg/mrcommons/man/calcProduction.html),
[`calcPriceAgriculture`](https://rdrr.io/pkg/mrcommons/man/calcPriceAgriculture.html)

## Author

David Chen

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("ValueProduction")
} # }
```
