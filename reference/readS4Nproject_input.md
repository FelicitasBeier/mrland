# readS4Nproject_input

reads in total bioenergy (1st gen, 2nd gen and residues) demand and co2
prices from IMAGE model for Sim4Nexus project

## Usage

``` r
readS4Nproject_input(subtype = "co2prices")
```

## Arguments

- subtype:

  IMAGE input to be read in: co2prices or bioenergy

## Value

magpie object at country-level resolution

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html)

## Author

Felicitas Beier

## Examples

``` r
if (FALSE) { # \dontrun{
a <- readSource("S4Nproject_input", convert = "onlycorrect", aggregate = FALSE)
} # }
```
