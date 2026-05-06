# Calculate global under-/overproduction

Calculate the difference between the global production and global
domestic_supply (corrective balance flow). The difference is the result
of imports not equaling exports, and because storage is not considered.
The calculated DomesticBalanceflow assures that production matches
domestic_supply. The goods come from nowhere and go to nowhere. The
numbers are usually decreased linearly and become zero in 2050.

## Usage

``` r
calcTradeBalanceflow()
```

## Value

global domestic balanceflows as MAgPIE object

## See also

[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html),
[`calcFAOmassbalance`](https://rdrr.io/pkg/mrcommons/man/calcFAOmassbalance.html)

## Author

Ulrich Kreidenweis, Xiaoxi Wang

## Examples

``` r
if (FALSE) { # \dontrun{
a <- calcTradeBalanceflow()
} # }
```
