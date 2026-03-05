#' @title calcIrrigUnitCost
#'
#' @description This function returns unit costs for irrigation expansion (new_totalCost)
#'              and maintenance or rehabilitation (rehab_hardwareCost)
#'
#' @return MAgPIE object
#'
#' @author Felicitas Beier
#'
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link{readInocencio2007}},
#' \code{\link{convertInocencio2007}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcIrrigUnitCost")
#' }
#'

calcIrrigUnitCost <- function() {

  # Read unit cost for irrigation by Inocencio et al. 2007
  data <- readSource("Inocencio2007")[, , c("new_totalCost", "rehab_hardwareCost")]
  getNames(data) <- c("expansion", "maintenance")

  # conversion: $2000 to $2017
  data <- toolConvertGDP(data,
                         unit_in = "constant 2000 US$MER",
                         unit_out = "constant 2017 US$MER",
                         replace_NAs = "no_conversion")

  return(list(x = data,
              weight = NULL,
              unit = "US$MER2017",
              description = "unit cost for irrigation expansion and maintenance"))
}
