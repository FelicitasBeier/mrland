#' @title calcGDPdeflator
#' @description calculates a iso-level deflator,
#' this is needed to run food demand and livestock regressions consistently
#' @param yearFrom year in "y2005" format
#' @param yearTo year in "y2005" format
#' @param currency "PPP" or "MER"
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author David Chen
#' @importFrom GDPuc toolConvertGDP

calcGDPdeflator <- function(yearFrom = 2017, yearTo = 2005, currency = "PPP") {

  reg <- getItems(calcOutput("GDP", scenario = "SSP2", unit = "constant 2005 Int$PPP", aggregate = FALSE), dim = 1)

  defl2017 <- new.magpie(cells_and_regions = reg,
                         years = NULL,
                         names = NULL,
                         fill = 1)

  if (currency == "PPP") {
    unit <- "Int$PPP"
  } else if (currency == "MER") {
    unit <- "US$MER"
  }

  defl2017 <- toolConvertGDP(defl2017,
                             unit_in = paste("constant", yearFrom, unit, sep = " "),
                             unit_out = paste("constant", yearTo, unit, sep = " "),
                             replace_NAs = c("linear", "no_conversion"))
  getNames(defl2017) <- NULL

  return(list(x = defl2017,
              weight = NULL,
              unit = "$/$",
              description = "GDP base year deflator"))

}
