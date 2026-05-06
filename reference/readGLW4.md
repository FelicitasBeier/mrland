# readGLW4

reads in Gridded Livestock of the World v4, downloaded from:
https://dataverse.harvard.edu/dataverse/glw_4

## Usage

``` r
readGLW4(subtype = "Da_Ct")
```

## Arguments

- subtype:

  Weighting method and livestock type:

  - Da: Dasymetric weighting informed by Random Forest

  - Aw: Areal weighting (distributed uniformly in each census)

    - `Ch`: Chicken

    - `Ct`: Cattle

    - `Pg`: Pigs

    - `Sh`: Sheep

    - `Gt`: Goats

    - `Ho`: Horse

    - `Dk`: ducks

    - `Bf`: Buffaloes

## Value

A gridded magpie object with gridded livstock of the world

## Author

David M Chen
