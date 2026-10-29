# readGFWLossByDriver

Reads the Global Forest Watch extract downloaded by
[`downloadGFWLossByDriver`](downloadGFWLossByDriver.md): tree cover loss
by country, year and driver, in Mha.

## Usage

``` r
readGFWLossByDriver()
```

## Value

magpie object, ISO country x year (2001-2025) x driver class, Mha

## See also

[`downloadGFWLossByDriver`](downloadGFWLossByDriver.md),
[`convertGFWLossByDriver`](convertGFWLossByDriver.md)

## Author

Michael Crawford

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("GFWLossByDriver")
} # }
```
