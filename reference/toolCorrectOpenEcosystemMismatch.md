# toolCorrectOpenEcosystemMismatch

Corrects mismatches in the classification of open ecosystems into
pasture and other between land-use initialisation (LUH) and ESA CCI.
Excess "other" land relative to LUH is shifted to "past", but only up to
the available LUH pasture surplus, to avoid spurious pasture expansion.

## Usage

``` r
toolCorrectOpenEcosystemMismatch(x, luIni)
```

## Arguments

- x:

  magpie object with "other" and "past" land-use classes. Cells and
  years already aligned with `luIni`

- luIni:

  magpie object from `calcOutput("LanduseInitialisation")` with
  aggregated "past" and "other" classes. Cells and years already aligned
  with `x`

## Value

`x` with the other/pasture mismatch corrected

## Author

Patrick v. Jeetze
