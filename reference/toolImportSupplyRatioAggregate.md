# toolImportSupplyRatioAggregate This function aggregates bilateral import supply ratios (amount exported from country X divided by domestic supply of country Y) from country to regional level, using two different weights for bilateral spatial dimensions 1.1 and 1.2, where dimension 1.1 is the exporters (no weight, simple summation), while dimension 1.2 are the importers, weighted by domestic supply. THIS ORDER IS IMPORTANT.

toolImportSupplyRatioAggregate This function aggregates bilateral import
supply ratios (amount exported from country X divided by domestic supply
of country Y) from country to regional level, using two different
weights for bilateral spatial dimensions 1.1 and 1.2, where dimension
1.1 is the exporters (no weight, simple summation), while dimension 1.2
are the importers, weighted by domestic supply. THIS ORDER IS IMPORTANT.

## Usage

``` r
toolImportSupplyRatioAggregate(x, rel, weight)
```

## Arguments

- x:

  magclass object that should be aggregated

- rel:

  relation matrix containing a region mapping.

- weight:

  weight to use for dimension 1.2

## Value

return: returns region aggregated import supply ratios

## Author

David M Chen
