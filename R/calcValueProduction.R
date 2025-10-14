#' @title calcValueProduction
#' @description calculates production value based on production and prices, only works for FAO dataset currently
#' @param cellular cellular or iso country values
#' @param datasource Options of the source of the price data: only FAO has country level data
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author David Chen
#' @seealso
#' \code{\link[mrcommons]{calcProduction}},
#' \code{\link[mrcommons]{calcPriceAgriculture}}
#' @examples
#' \dontrun{
#' calcOutput("ValueProduction")
#' }
#'
calcValueProduction <- function(datasource = "FAO", cellular = TRUE) {

if (datasource == "FAO") {

  if (cellular == TRUE) {
  prod <- calcOutput("Production", products = "kcr",
                     cellular = TRUE,
                     aggregate = FALSE)[, , "dm"]
  } else {
  prod <- calcOutput("Production", products = "kcr",
                     cellular = FALSE, aggregate = FALSE)[, , "dm"]
  }

prod <- collapseNames(prod)

prices <- calcOutput("PriceAgriculture", datasource = datasource, aggregate = FALSE)
prices <- collapseNames(prices)


common  <- intersect(getNames(prices), getNames(prod))
regions <- intersect(getItems(prices, dim = 1.1), getItems(prod, dim = 1.1))

x <- prod[regions, , common] * prices[regions, , common]
} else {
  stop("Only FAO datasource has country level prices")
}
return(list(x = x,
            weight = NULL,
            unit = "US$2017/tDM",
            description = "Crop Production Value",
            min = 0,
            max = Inf,
            isocountries = !cellular))
}
