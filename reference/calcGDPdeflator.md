# calcGDPdeflator

calculates a iso-level deflator, this is needed to run food demand and
livestock regressions consistently

## Usage

``` r
calcGDPdeflator(yearFrom = 2017, yearTo = 2005, currency = "PPP")
```

## Arguments

- yearFrom:

  year in "y2005" format

- yearTo:

  year in "y2005" format

- currency:

  "PPP" or "MER"

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## Author

David Chen
