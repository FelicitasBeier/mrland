# convertGFWLossByDriver

Moves the GFW data from the GADM country set to madrat's ISO country
set: Kosovo, Northern Cyprus and the Akrotiri and Dhekelia base areas
are added to Serbia and Cyprus, GADM's placeholder codes for disputed
territories are dropped and countries GFW does not list are filled with
zero.

## Usage

``` r
convertGFWLossByDriver(x)
```

## Arguments

- x:

  magpie object as returned by
  [`readGFWLossByDriver`](readGFWLossByDriver.md)

## Value

magpie object on madrat's ISO country set, Mha

## See also

[`readGFWLossByDriver`](readGFWLossByDriver.md)

## Author

Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("GFWLossByDriver", convert = TRUE)
} # }
```
