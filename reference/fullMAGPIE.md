# fullMAgPIE

Function that produces the regional data set for running the MAgPIE
model.

## Usage

``` r
fullMAGPIE(rev = numeric_version("0.1"), dev = "")
```

## Arguments

- rev:

  data revision which should be used as input (numeric_version).

- dev:

  For developing purposes, apply changes as per dev flag

## See also

[`readSource`](https://rdrr.io/pkg/madrat/man/readSource.html),
[`getCalculations`](https://rdrr.io/pkg/madrat/man/getCalculations.html),
[`calcOutput`](https://rdrr.io/pkg/madrat/man/calcOutput.html)

## Author

Jan Philipp Dietrich, Benjamin Leon Bodirsky, Florian Humpenoeder, Edna
J. Molina Bacca

## Examples

``` r
if (FALSE) { # \dontrun{
retrieveData("MAGPIE", rev = numeric_version("12"),
             mainfolder = "pathtowhereallfilesarestored")
} # }
```
