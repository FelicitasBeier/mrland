#' @title calcBMIshr
#'
#' @description estimates population based on BMI share
#'
#' @param convert Use raw data or interpolated data. Raw data should only be used for regressions.
#'
#' @return List with a magpie object
#' @author Benjamin Leon Bodirsky
#' @seealso
#' \code{\link{readNCDrisc}},
#' \code{\link{calcIntake}}
#'
#' @examples
#' \dontrun{
#' calcOutput("BMIshr")
#' }
#'
#' @importFrom magclass getRegions dimSums

calcBMIshr <- function(convert = TRUE) {
  ### Adult

  x <- readSource("NCDrisc", subtype = "BMI_shr", convert = FALSE)
  mapping <- toolGetMapping(type = "sectoral", name = "NCDriscBMIshr2Lutz.csv", where = "mappingfolder")
  x <- toolAggregate(x, rel = mapping, from = "NCDrisc", to = "lutz", dim = 3.1)
  mapping <- toolGetMapping(type = "sectoral", name = "BMIgroup_adultBMI.csv", where = "mappingfolder")
  adult <- toolAggregate(x, rel = mapping, from = "adultBMI", to = "BMIgroup", dim = 3.3, weight = NULL, partrel = TRUE)

  ### underaged

  x <- readSource("NCDrisc", subtype = "BMI_shr_underaged", convert = FALSE)

  ### aggregate to age groups, use age 5 for 0--4
  relevant <- c("age5", "age6", "age7", "age8", "age9", "age10", "age11", "age12", "age13", "age14")

  x <- x[, , relevant]
  weight <- x * 0 + 1  ### assume equal weighting within age groups
  mapping <- toolGetMapping(type = "sectoral", name = "NCDrisc2Lutz.csv", where = "mappingfolder")
  x <- toolAggregate(x, rel = mapping, from = "NCDrisc", to = "lutz", dim = 3.1, weight = weight, partrel = TRUE)
  mapping <- toolGetMapping(type = "sectoral", name = "BMIgroup_underagedBMI.csv", where = "mappingfolder")
  x <- toolAggregate(x, rel = mapping, from = "underagedBMI", to = "BMIgroup", dim = 3.3, weight = NULL, partrel = TRUE)
  underaged <- add_columns(x, addnm = c("mediumhigh"), dim = 3.3)
  underaged[, , c("mediumhigh")] <- 0
  ###

  out <- mbind(adult, underaged)
  out <- out[, , getNames(adult, dim = 3)] ### right order

  if (convert) {
    bmi <- out
    withdata <- getItems(bmi, 1.1)
    bmi2 <- toolCountryFill(bmi, fill = NA)
    bmi2 <- add_columns(bmi2, dim = 2.1, addnm = c("y1965", "y1970"))
    bmi2 <- bmi2[, sort(getYears(bmi2)), ]

    regression <- readSource("Bodirsky2018", convert = FALSE)

    gdpPC <- collapseNames(calcOutput("GDPpc", scenario = "SSP2", aggregate = FALSE)[, , ])

    bmiRegr <- collapseNames(regression[, , "intercept"]
                             + regression[, , "saturation"] * gdpPC / (regression[, , "halfsaturation"] + gdpPC))
    bmiRegr <- time_interpolate(bmiRegr, interpolated_year = getYears(bmi2), integrate_interpolated_years = FALSE)

    mapping <- toolGetMapping(type = "sectoral", name = "Lutz2agegroups.csv", where = "mappingfolder")


    bmi2 <- bmi2 * NA

    for (agegroup in c("underaged", "working", "retired")) {
      ages <- mapping[mapping$agegroups == agegroup, 1]
      bmi2[, , "verylow"][, , ages] <- bmiRegr[, , agegroup][, , "low"] * bmiRegr[, , agegroup][, , "lowsplit"]
      bmi2[, , "low"][, , ages] <- bmiRegr[, , agegroup][, , "low"] * (1 - bmiRegr[, , agegroup][, , "lowsplit"])
      bmi2[, , "medium"][, , ages] <- ((1 - bmiRegr[, , agegroup][, , "low"] - bmiRegr[, , agegroup][, , "high"])
                                       * (1 - bmiRegr[, , agegroup][, , "mediumsplit"]))
      bmi2[, , "mediumhigh"][, , ages] <- ((1 - bmiRegr[, , agegroup][, , "low"] - bmiRegr[, , agegroup][, , "high"])
                                           * (bmiRegr[, , agegroup][, , "mediumsplit"]))
      bmi2[, , "high"][, , ages] <- bmiRegr[, , agegroup][, , "high"] * (1 - bmiRegr[, , agegroup][, , "highsplit"])
      bmi2[, , "veryhigh"][, , ages] <- bmiRegr[, , agegroup][, , "high"] * bmiRegr[, , agegroup][, , "highsplit"]
    }

    calib <- bmi[, "y1975", ] - bmi2[withdata, "y1975", ]

    bmi2[withdata, getYears(bmi), ] <- bmi[withdata, getYears(bmi), ]

    bmi2[withdata, c("y1965", "y1970"), ] <- bmi2[withdata, c("y1965", "y1970"), ] + setYears(calib, NULL)
    # in case that calibration created negative values or values above one,
    # remove them and add them to the middle category
    bmi2[bmi2 < 0] <- 0.000001
    bmi2[bmi2 > 1] <- 1
    bmi2[, , "medium"] <- bmi2[, , "medium"] + (1 - dimSums(bmi2, dim = 3.3))

    out <- bmi2

  }

  weight <- collapseNames(calcOutput("Demography", aggregate = FALSE, education = FALSE)[, , "SSP2"])
  weight <- time_interpolate(weight, interpolated_year = getYears(out), extrapolation_type = "constant")

  return(list(x = out,
              weight = weight,
              unit = "capita/capita",
              description = "Share of population belonging to a BMI group.",
              isocountries = convert))
}
