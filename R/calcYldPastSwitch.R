#' calcYldPastSwitch
#'
#' Dummy file for regional pasture yield switch (f14_yld_past_switch)
#'
#' @return Dummy file for regional pasture yield switch
#' @author Alexandre Köberle, Patrick Rein
#' @seealso \code{\link[madrat]{readSource}}, \code{\link[madrat]{calcOutput}}
#' @importFrom magclass new.magpie
#' @importFrom madrat toolGetMapping

calcYldPastSwitch <- function() {
  isoCountry  <- toolGetMapping("iso_country.csv", where = "mrland")
  isoCountry1 <- as.vector(isoCountry[, "x"])
  names(isoCountry1) <- isoCountry[, "X"]
  x <- new.magpie(cells_and_regions = isoCountry1, years = seq(1965, 2150, by = 5),
                  names = NULL, fill = 0.25)

  return(list(x = x,
              weight = x,
              unit = "-",
              description = "Dummy file for regional pasture yield switch"))
}
