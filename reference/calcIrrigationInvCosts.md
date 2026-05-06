# calcIrrigationInvCosts

This function calculates irrigation investment costs for each country
until the year 2050. Values linearly converge towards the value of
Germany (1995) by 2050.

## Usage

``` r
calcIrrigationInvCosts()
```

## Value

MAgPIE object

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`readWBirrigation`](readWBirrigation.md),
[`convertWBirrigation`](convertWBirrigation.md)

## Author

Nele Steinmetz, Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("IrrigationInvCosts")
} # }
```
