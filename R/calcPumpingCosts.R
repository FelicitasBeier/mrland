#' @title calcPumpingCosts
#' @description provides costs of pumping irrigation water
#' @return A magpie object at iso level for all years with information on pumping costs
#' @author Vartika Singh
#' #' @seealso \code{\link[madrat]{readSource}}, \code{\link[madrat]{calcOutput}}
#' @importFrom magclass new.magpie
#' @importFrom utils read.csv2
#' @examples
#' \dontrun{
#' calcOutput("PumpingCosts")
#' }
#'
calcPumpingCosts <- function() {

  isoCountry  <- toolGetMapping("iso_country.csv", where = "mrland")
  isoCountry1 <- as.vector(isoCountry[, "x"])
  names(isoCountry1) <- isoCountry[, "X"]
  x <- new.magpie(cells_and_regions = isoCountry1, years = seq(1995, 2100, by = 5), names = NULL, fill = 0)

  # Assigning a value of 0.005 cents for India
  x["IND", , ] <- 0.005

  weight <- x
  weight[, , ] <- 1

  return(list(x            = x,
              weight       = weight,
              unit         = "USD per million cubic meters",
              description  = "costs of pumping irrigation water from Cornish et.al., 2004"))

}
