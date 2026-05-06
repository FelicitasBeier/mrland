# readREMIND

Reads in a reporting mif file from REMIND

## Usage

``` r
readREMIND(subtype)
```

## Arguments

- subtype:

  A string composed of three items: unit, revision and indicator. Unit
  can be either "intensive" or "extensive", revision is the input data
  revision, and indicator is the name of the REMIND indicator

## Value

MAgPIE object with regional aggregation of REMIND H12

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

David Klein

## Examples

``` r

if (FALSE) { # \dontrun{
readSource("REMIND",aggregate=FALSE)
} # }
```
