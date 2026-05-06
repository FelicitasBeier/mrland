# calcTradeStdDevHelper

Calculate rolling standard-deviation statistics for trade ratios
observed over a window of years.

## Usage

``` r
calcTradeStdDevHelper(
  dataIn,
  sdYears,
  groupVars = c("ex", "im", "ItemCodeItem"),
  yearCol = "Year",
  valueCol = ".value"
)
```

## Arguments

- dataIn:

  data.frame or tibble; tidy table with columns specified in
  \`group_vars\`, \`year_col\` and \`value_col\`.

- sdYears:

  integer; window length (number of years) for the rolling standard
  deviation (must be \>= 1).

- groupVars:

  character; grouping columns. Default: c("ex", "im", "ItemCodeItem").

- yearCol:

  character; name of the year column. Default: "Year".

- valueCol:

  character; name of the numeric value column. Default: ".value".

## Value

tibble with one row per group and columns: meansd\<years\>,
maxsd\<years\>, minsd\<years\>.

## Details

Given a tidy data.frame/tibble with columns for exporter (ex), importer
(im), item (ItemCodeItem), year and a numeric value column (default
\`.value\`), this function computes a rolling standard deviation over
\`years\` for each group and returns the mean, max and min of that
rolling standard deviation per group.

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
ratiodf <- as.data.frame(collapseNames(ratio), rev = 2)
toolTradeStdDev(ratiodf, years = 5)
} # }
```
