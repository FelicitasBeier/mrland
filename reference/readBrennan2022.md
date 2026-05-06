# readBrennan2022

Reads Critical Connectivity Areas as defined in Brennan, A., Naidoo, R.,
Greenstreet, L., Mehrabi, Z., Ramankutty, N., & Kremen, C. (2022).
Functional connectivity of the world’s protected areas. Science,
376(6597), 1101–1104. https://doi.org/10.1126/science.abl8974 Protected
areas (2020) and Key Biodiversity Areas/Global Safet Net areas were
masked at a spatial resolution of 10 arc seconds before aggregating the
data to 0.5°.

## Usage

``` r
readBrennan2022(subtype = "KBA_GSN_masked")
```

## Arguments

- subtype:

  Defines whether land area covered by Critical Connectivity Areas has
  been masked by other conservation priority data. If Key Biodiversity
  Areas have only been masked the option is `"KBA_masked"`. With
  `"KBA_GSN_masked"`, land area covered by the Global Safety Net
  (distinct species assemblages cluster) is also masked. This is useful
  for complementary scenario building.

## Value

Returns magpie objects with the land area covered by Critical
Connectivity areas that is NOT already covered by Key Biodiversity Areas
or the Global Safety Net (distinct species assemblages cluster) and was
unprotected in 2020.

## Author

Patrick v. Jeetze

## Examples

``` r
if (FALSE) { # \dontrun{
readSource("Brennan2022", convert = "onlycorrect")
} # }
```
