#' @title calcTradeBilateralScenarioAdj
#'
#' @description Creates dummy trade scenario adjustments file filld with 0's,
#' actual changes currently implemented in the model on the import supply ratio
#' (imports/domestic supply)
#' @param magYears logical or character vector. If TRUE, filters to t_all 5-year
#' time steps using findset("t_all").
#' @return dummy file with 0's to befilled in MAgPIE
#' @author David M Chen
#'
calcTradeBilateralScenarioAdj <- function(magYears = TRUE) {

  tm <- calcOutput("TradeBilateralFAOHarmonized", aggregate = FALSE, yearly = TRUE)
  tm[] <- 0

  if (magYears) {
    tAll <- findset("t_all")
    tm <- tm[, intersect(getYears(tm), tAll), ]
  }
  tm <- toolHoldConstantBeyondEnd(tm)
  #replace set names order from im.ex to ex.im
  getSets(tm)[c(1, 2)] <- getSets(tm)[c(2, 1)]

  return(list(x = tm,
              weight = NULL,
              unit = "percentage",
              description = "changes to countries import supply ratio"))
}
