# calcBiocharProd

calculates biochar production

## Usage

``` r
calcBiocharProd(
  datasource,
  rev = numeric_version("0.1"),
  biocharsystems = "disagg",
  unit = "EJ/yr"
)
```

## Arguments

- datasource:

  source to be used

- rev:

  data revision the output will be produced for (numeric_version).

- biocharsystems:

  selects the biochar production systems for which data is returned.
  Options are disagg which returns data for all available biochar
  production systems separately, and agg which returns the total of all
  biochar production systems.

- unit:

  (available: EJ/yr and Mt/yr)

## Value

magpie object with results on country level, weight on country level,
unit and description.

## See also

[`calc2ndBioDem`](calc2ndBioDem.md), [`calcGHGPrices`](calcGHGPrices.md)

## Author

Isabelle Weindl

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("BiocharProd")
} # }
```
