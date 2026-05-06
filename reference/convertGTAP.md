# convertGTAP

Converts GTAP data to fit to the common country list. Weighting is done
by using the Imports and Exports from FAO. NOW NEW WEIGHTING

## Usage

``` r
convertGTAP(x, subtype)
```

## Arguments

- x:

  MAgPIE object contains GTAP data

- subtype:

  The GTAP subtype: VIWS, VIMS VXWD, VXMD, VOA, VOM

## Value

Converted GTAP Data

## Author

Xiaoxi Wang

## Examples

``` r
if (FALSE) { # \dontrun{
x <- ReadSource("GTAP", "GTAP7_VIMS")
} # }
```
