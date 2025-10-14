#' @title calcUrbanLandFuture
#' @description Urban land in Mha on 0.5deg grid
#' @param cellular TRUE for results on 0.5 degree grid.
#' @param timestep 5year or yearly
#' @param cells    (deprecated) lpjcell (67420 cells)
#' @param subtype  where the data source comes from ("LUH3" or "Gao")
#' @return List of magpie objects with results on 0.5deg grid level, weights NULL, unit and description.
#' @author David Chen, Patrick v. Jeetze, Felicitas Beier
#' @importFrom magpiesets findset
#' @importFrom mstools toolHoldConstant
#' @importFrom magclass nregions setCells getCells

calcUrbanLandFuture <- function(timestep = "5year", subtype = "LUH3",
                                cells = "lpjcell", cellular = TRUE) {

  if (subtype == "LUH3") {
    past <- calcOutput("LanduseInitialisation",
      cellular = TRUE, nclasses = "seven",
      selectyears = "past_til2020", input_magpie = FALSE, aggregate = FALSE
    )
    past <- past[, c("y1995", "y2000", "y2005", "y2010", "y2015", "y2020"), "urban"]
    past <- add_columns(past, addnm = c("SSP2", "SSP3", "SSP4", "SSP5"), dim = 3.1)
    past[, , 2:5] <- past[, , 1]
    getNames(past)[1] <- "SSP1"

    out <- readSource("LUH2UrbanFuture", convert = "onlycorrect")
    out <- out[, c(2015:2020), invert = TRUE]

    if (timestep == "5year") {
      out <- out[, paste0("y", seq(2025, 2100, 5)), ]
      out <- toolHoldConstant(out, paste0("y", seq(2105, 2150, 5)))
      out <- mbind(setCells(past, getCells(out)), out)
      names(dimnames(out)) <- c("x.y.iso", "t", "data")
    } else if (timestep == "yearly") {
      out <- toolHoldConstant(out, paste0("y", c(2101:2150)))
      past <- time_interpolate(past, interpolated_year = 1995:2020)
      out <- mbind(setCells(past, getCells(out)), out)
      names(dimnames(out)) <- c("x.y.iso", "t", "data")
    }

  } else if (subtype == "Gao") {

    out <- readSource("UrbanLandGao", convert = FALSE)

    if (timestep == "5year") {
      years <- seq(1995, 2100, 5)
      out <- time_interpolate(out, interpolated_year = years, integrate_interpolated_years = TRUE)
      out <- time_interpolate(out,
        extrapolation_type = "constant",
        interpolated_year = seq(2105, 2150, 5), integrate_interpolated_years = TRUE
      )
      # set any interpolated values to 0
      out[out < 0] <- 0
    }

  } else {
    stop("Not a Valid Subtype")
  }

  # make all years
  yHarm <- seq(1995, 2025, 5)
  for (i in c("SSP1", "SSP3", "SSP4", "SSP5")) {
    out[, yHarm, i] <- out[, yHarm, "SSP2"]
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    isocountries = (!cellular & (nregions(out) != 1)),
    description = "Amount of Urban land expansion for various SSPs"
  ))
}
