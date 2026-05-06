# calcTradeMargin

calculate total value of trade margins from GTAP dataset

## Usage

``` r
calcTradeMargin(
  gtap_version = "GTAP9",
  bilateral = FALSE,
  producer_price = "FAOini"
)
```

## Arguments

- gtap_version:

  type of GTAP data version

  - `GTAP7`

  - `GTAP8`

  - `GTAP9`

- bilateral:

  whether bilateral trade margin should be calculated

- producer_price:

  which producer price should be used

## Value

Trade margins as an MAgPIE object

## Author

Xiaoxi Wang

## Examples

``` r
if (FALSE) { # \dontrun{
x <- calcTradeMargin("GTAP7")
} # }
```
