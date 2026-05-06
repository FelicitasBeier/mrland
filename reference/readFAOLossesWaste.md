# Read in data on food losses and waste from FAO for several commodity groups

Data from Annex 4 of the following FAO study: FAO. 2011. Global food
losses and food waste – Extent, causes and prevention. Rome
(http://www.fao.org/3/a-i2697e.pdf)

## Usage

``` r
readFAOLossesWaste(subtype)
```

## Arguments

- subtype:

  Steps of the food supply chain where food losses and waste occur.
  Available types are:

  - `Consumption`: consumption level

## Value

magpie object of food waste percentages for several commodity groups

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl

## Examples

``` r

if (FALSE)  a <- readSource(type="FAOLossesWaste",subtype="Consumption")
 # \dontrun{}
```
