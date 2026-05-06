# calcGHGPrices

reads in GHG prices from past runs

## Usage

``` r
calcGHGPrices(
  emissions = "pollutants",
  datasource = "REMMAG",
  rev = numeric_version("0.1")
)
```

## Arguments

- emissions:

  which type of emissions shall be returned. ghg just returns n2o, ch4
  and co2, pollutants a longer list including also air pollutants

- datasource:

  REMIND for prices from R2M4 coupled runs, REMMAG for old coupled runs,
  SSPResults for prices from the SSP scenarios from the IIASA database,
  SSP_and_REM for a combination of REMIND and SSPResults

- rev:

  data revision the output will be produced for (numeric_version).

## Value

list of magpie object with results on country level, weight on country
level, unit and description.

## See also

[`readSSPResults`](https://rdrr.io/pkg/mrcommons/man/readSSPResults.html)

## Author

David Chen, Benjamin Leon Bodirsky, David Klein

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("GHGPrices")
} # }
```
