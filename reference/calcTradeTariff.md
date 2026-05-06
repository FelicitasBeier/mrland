# calcTradeTariff

calculate tarde tariffs from GTAP dataset

## Usage

``` r
calcTradeTariff(
  gtap_version = "GTAP9",
  type_tariff = "total",
  bilateral = FALSE
)
```

## Arguments

- gtap_version:

  type of GTAP data version

  - `GTAP7`

  - `GTAP8`

- type_tariff:

  which producer price should be used

  - `type_tariff`

- bilateral:

  calculates whether tariffs should be bilateral

## Value

Trade tariffs as an MAgPIE object

## Author

Xiaoxi Wang, David M Chen

## Examples

``` r
    if (FALSE) { # \dontrun{
    x <- calcTradeTariff("GTAP7")
    } # }
```
