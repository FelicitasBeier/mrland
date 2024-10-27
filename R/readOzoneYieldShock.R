#' @title readOzoneYieldShock
#' @description read Ozone Yield Shock
#' Data from the EAT-Lancet deepdive on Ozone shock effects on crop yields.
#' @return MAgPIE object with country level yield shock data for year 2050.
#' @author Jake Tommey
#' @examples
#' \dontrun{
#' readSource("OzoneShock", convert = "onlycorrect")
#' }
#' @importFrom readxl read_xlsx

readOzoneYieldShock <- function() {
  # read-in file for ozone shocks
  x <- readxl::read_excel("Ag_CCShocks.xlsx", sheet = "Ozone")
  # convert to magclass object
  x <- as.magpie(x[, 1:4], spatial = 1)
  getYears(x) <- 2050
  return(x)
}
