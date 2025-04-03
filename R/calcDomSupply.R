#' Calculate regional domestic supply for trade module
#'
#' Simple function to get the domestic supply for the trade module
#' @title calcTradeMargin
#' @return magpie object domestic supply
#' @param nutrient "dm"
#' @param extend extend into magpie future years
#' @author David Chen
#' @examples
#' \dontrun{
#' a <- calcDomSupply()
#' }
#'
calcDomSupply <- function(nutrient = "dm", extend = TRUE) {

  out <- collapseNames(calcOutput("FAOmassbalance", aggregate = FALSE)[, , nutrient][, , "domestic_supply"])

  if (extend) {
    out <-  toolHoldConstantBeyondEnd(out)
  }

  return(list(x = out,
              weight = NULL,
              unit = "Mt",
              description = "Domestic supply"))
}
