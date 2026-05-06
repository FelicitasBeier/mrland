# calcEFch4Rice

emission factors for methane from flooded rice fields, depending on
phyiscal area or area harvested. The emission factors were calculated
based on FAOSTAT estimates due to lack of all necessary parameters in
the IPCC Guidelines

## Usage

``` r
calcEFch4Rice(physical = TRUE)
```

## Arguments

- physical:

  if true physical area, if false area harvested

## Value

List of magpie objects with results on country level, weight on country
level, unit and description.

## See also

[`calcEFch4AWMS`](calcEFch4AWMS.md),
[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html)

## Author

Benjamin Leon Bodirsky

## Examples

``` r
if (FALSE) { # \dontrun{
calcOutput("EFch4Rice")
} # }
```
