# calcTradeBilateralScenarioAdj

Creates dummy trade scenario adjustments file filld with 0's, actual
changes currently implemented in the model on the import supply ratio
(imports/domestic supply)

## Usage

``` r
calcTradeBilateralScenarioAdj(magYears = TRUE)
```

## Arguments

- magYears:

  logical or character vector. If TRUE, filters to t_all 5-year time
  steps using findset("t_all").

## Value

dummy file with 0's to befilled in MAgPIE

## Author

David M Chen
