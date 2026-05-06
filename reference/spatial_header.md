# Tool: spatial_header

Given a regionmapping (mapping between ISO countries and regions) the
function calculates a 0.5 degree spatial header for 0.5 degree magclass
objects

## Usage

``` r
spatial_header(mapping)
```

## Arguments

- mapping:

  Either a path to a mapping or an already read-in mapping as
  data.frame.

## Value

A vector with 59199 elements

## See also

[`regionscode`](https://rdrr.io/pkg/madrat/man/regionscode.html)

## Author

Jan Philipp Dietrich

## Examples

``` r
if (FALSE) { # \dontrun{
spatial_header("regionmappingMAgPIE.csv")
} # }
```
