# Read in data from the NIN recommendations

Read in data from the NIN recommendations

## Usage

``` r
readNIN(subtype)
```

## Arguments

- subtype:

  Type of NIN data that should be read. Available types are:

  - `cons_data`: Consumption analysis ("NIN_cons_data.csv")

## Value

magpie object containing NIN data

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Isabelle Weindl, Jan Philipp Dietrich

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource(type = "NIN", subtype = "cons_data")
} # }
```
