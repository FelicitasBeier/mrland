# calcEmisNitrogenPast

Emission factors from cropland soils.

## Usage

``` r
calcEfNSoil(method = "IPCC_reg")
```

## Arguments

- method:

  If IPCC, using the ipcc emission factors as share of applied N inputs.
  If Nloss, as share of cropland budget surplus.

## Value

List of magpie object with results on country level, weight on country
level, unit and description.

## Author

Benjamin Leon Bodirsky

## Examples

``` r

if (FALSE) { # \dontrun{
calcOutput("EmisNitrogenPast")
} # }
```
