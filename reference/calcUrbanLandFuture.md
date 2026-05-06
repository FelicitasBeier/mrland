# calcUrbanLandFuture

Urban land in Mha on 0.5deg grid

## Usage

``` r
calcUrbanLandFuture(
  timestep = "5year",
  subtype = "LUH3",
  cells = "lpjcell",
  cellular = TRUE
)
```

## Arguments

- timestep:

  5year or yearly

- subtype:

  where the data source comes from ("LUH3" or "Gao")

- cells:

  (deprecated) lpjcell (67420 cells)

- cellular:

  TRUE for results on 0.5 degree grid.

## Value

List of magpie objects with results on 0.5deg grid level, weights NULL,
unit and description.

## Author

David Chen, Patrick v. Jeetze, Felicitas Beier
