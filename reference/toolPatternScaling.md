# toolPatternScaling

This tool scales time series based on the approach used in the
magpiemodel yield module.

## Usage

``` r
toolPatternScaling(
  scen,
  scenMean,
  refMean,
  refYear = "y2010",
  variation = "yieldCalibMAG"
)
```

## Arguments

- scen:

  time series of the scenario

- scenMean:

  mean of scenario time series

- refMean:

  mean of reference time series

- refYear:

  Reference year

- variation:

  'yieldCalibMAG' (default); to be implemented: 'jensPaper'

## Value

scaled data in magclass format

## Author

Kristine Karstens
