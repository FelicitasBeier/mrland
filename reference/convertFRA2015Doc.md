# Convert FRA2015Doc data

Convert FRA2015Doc data

## Usage

``` r
convertFRA2015Doc(x, subtype)
```

## Arguments

- x:

  MAgPIE object containing original values coming from read function

- subtype:

  The data table type, e.g.: forest_area

## Value

Data as MAgPIE object

## See also

[`readFRA2015Doc`](readFRA2015Doc.md),
[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html),

## Author

Abhijeet Mishra

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("FRA2015Doc", "forest_area", convert = TRUE)
} # }
```
