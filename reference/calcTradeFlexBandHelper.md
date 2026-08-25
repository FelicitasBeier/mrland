# calcTradeFlexBandHelper

Calculate the rolling-range flexibility band of trade ratios observed
over a window of years.

## Usage

``` r
calcTradeFlexBandHelper(
  dataIn,
  windowYears,
  groupVars = c("ex", "im", "ItemCodeItem"),
  yearCol = "Year",
  valueCol = ".value"
)
```

## Arguments

- dataIn:

  data.frame or tibble; tidy table with columns specified in
  \`groupVars\`, \`yearCol\` and \`valueCol\`.

- windowYears:

  integer; window length (number of years) for the rolling range (must
  be \>= 1).

- groupVars:

  character; grouping columns. Default: c("ex", "im", "ItemCodeItem").

- yearCol:

  character; name of the year column. Default: "Year".

- valueCol:

  character; name of the numeric value column. Default: ".value".

## Value

tibble with one row per group and columns: mean\<years\>, max\<years\>,
min\<years\>.

## Details

Given a tidy data.frame/tibble with columns for exporter (ex), importer
(im), item (ItemCodeItem), year and a numeric value column (default
\`.value\`), this function computes the rolling RANGE (max - min) over
\`windowYears\` for each group and returns the mean, max and min of that
rolling range per group. The mean variant (mean over windows of the
rolling range) is used by the model as the flexibility window; it is
empirically monotone in window length on the historical data. The max
variant is provably monotone (range is a monotone set function) and
wider; min is the narrowest. All three are written to the input file.

## Author

David M Chen

## Examples

``` r
if (FALSE) { # \dontrun{
ratiodf <- as.data.frame(collapseNames(ratio), rev = 2)
calcTradeFlexBandHelper(ratiodf, windowYears = 5)
} # }
```
