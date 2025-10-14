#' @title calcIrrigationInvCosts
#'
#' @description This function calculates irrigation investment costs for each country until
#'              the year 2050. Values linearly converge towards the value of Germany (1995)
#'              by 2050.
#'
#' @return MAgPIE object
#'
#' @author Nele Steinmetz, Felicitas Beier
#'
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link{readWBirrigation}},
#' \code{\link{convertWBirrigation}}
#'
#' @examples
#' \dontrun{
#' calcOutput("IrrigationInvCosts")
#' }
#'
#' @import magclass
#'
calcIrrigationInvCosts <- function() {

  # Investment costs for expanding irrigation infrastructure in 1000 US$ per hectare
  # Note: in World Bank report (William I. Jones, 1991), p. 98, given in US$ per hectare
  # But: csv table has transformed values (in 1000 US$ per hectare)
  wBirrigation           <- readSource("WBirrigation")
  getYears(wBirrigation) <- NULL

  # irrigation cost constant until 2015
  data              <- new.magpie(cells_and_regions = getItems(wBirrigation, dim = 1),
                                  years = 1995:2050,
                                  names = "ad_unit_cost",
                                  fill = NA)
  data[, 1995:2050, ] <- wBirrigation

  # Transform fromm 1000 USD$ per hectare to US$ per hectare
  data <- data * 1000

  # conversion: $1995 to $2004
  data <- toolConvertGDP(data,
                         unit_in = "constant 1995 US$MER",
                         unit_out = "constant 2017 US$MER",
                         replace_NAs = "no_conversion")

  # from 2015 onwards, data converges to value of Germany until 2050
  dataDEU <- new.magpie(cells_and_regions = getItems(wBirrigation, dim = 1),
                        years = 1995:2050,
                        names = "ad_unit_cost",
                        fill = as.numeric(data["DEU", "y1995", ]))

  data <- convergence(origin = data, aim = dataDEU,
                      start_year = "y2015", end_year = "y2050",
                      type = "linear")

  # expand this value until 2150
  timeExtend <- paste0("y", seq(2055, 2150, 5))
  data       <- time_interpolate(data, timeExtend,
                                 extrapolation_type = "constant",
                                 integrate_interpolated_years = TRUE)
  getNames(data) <- NULL

  # aggregation weight and transform from 1000ha to Mha
  w <- readSource("FAO_online", subtype = "Land",
                  convert = TRUE)[, 1995, "6610|Agricultural land.Area_(1000_ha)"] * 1000

  return(list(x = data,
              weight = w,
              unit = "US$MER2017",
              description = "unit cost for irrigation expansion"))
}
