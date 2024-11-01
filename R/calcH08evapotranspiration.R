#' @title calcH08evapotranspiration
#' @description Calc evapotranspiration data for SSP cenarios in mm/month
#' @param subtype Switch between different inputs
#' @return magpie object in cellular resolution
#' @author Marcos Alves
#'
#' @examples
#' \dontrun{
#' calcOutput("H08evapotranspiration", subtype = "H08:mri-esm2-0")
#' }
#'
#' @import madrat
#' @import magclass
#' @importFrom magpiesets findset
#'

calcH08evapotranspiration <- function(subtype = "H08:mri-esm2-0") {

  time <- findset("time")
  timePast <- time[1:match("y2010", time)]
  timeFuture <- time[match("y2015", time):match("y2100", time)]

  y <-  readSource("H08evapotranspiration", subtype = paste0(subtype, ":", "historical"),
                   convert = FALSE)[, timePast, ]
  scenarios <- c("ssp126", "ssp370", "ssp585") # Current ISIMIP3bv2 scenarios
  x <- list()
  for (scenario in scenarios) {
    x[[scenario]] <-  mbind(setNames(y, scenario), readSource("H08evapotranspiration",
                                                              subtype = paste0(subtype, ":", scenario),
                                                              convert = FALSE)[, timeFuture, ])
  }
  x <- collapseNames(mbind(x))
  # from kg m-2 s-1 (https://protocol.isimip.org/#output-data -> evap)
  # to mm/day (https://gitlab.pik-potsdam.de/lpjml/LPJmL_internal/-/wikis/Input#climate)
  # (https://www.researchgate.net/post/How-do-I-convert-ERA-Interim-precipitation-estimates-from-kg-m2-s-to-mm-day)
  x <- x * 86400

  return(list(
    x = x,
    weight = NULL,
    unit = "mm/day",
    description = "Total evapotranspiration",
    isocountries = FALSE
  ))
}
