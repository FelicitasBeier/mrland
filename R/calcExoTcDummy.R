#' calcExoTcDummy
#'
#' Dummy file for regional exogenous tau path
#'
#' @return Dummy file for regional exogenous tau path
#' @author Florian Humpenoeder
#' @seealso \code{\link[madrat]{readSource}}, \code{\link[madrat]{calcOutput}}
#' @importFrom magclass new.magpie
#' @importFrom utils read.csv2

calcExoTcDummy <- function() {
  isoCountry  <- toolGetMapping("iso_country.csv", where = "mrland")
  isoCountry1 <- as.vector(isoCountry[, "x"])
  names(isoCountry1) <- isoCountry[, "X"]
  x <- new.magpie(cells_and_regions = isoCountry1, years = seq(1995, 2150, by = 5),
                  names = c("crop", "pastr"), fill = 0)

  return(list(x = x,
              weight = NULL,
              unit = "-",
              description = "Dummy file for regional exogenous tau path",
              note = "All values in the file are set to 0 if a new regional setup is used"))
}
