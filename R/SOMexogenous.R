#' @title calcSOMexogenous
#' @description
#' Uses an exogenous trajectory of Soil organic matter loss nitrogen release
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @author Benjamin Leon Bodirsky
#' @seealso \code{\link[madrat]{calcOutput}}
#' @examples
#' \dontrun{
#' calcOutput("SOMexogenous")
#' }
#' @export
#' @importFrom mstools toolHoldConstantBeyondEnd
#' @importFrom magpiesets findset

calcSOMexogenous <- function() {

  input <- calcOutput("SOMlossN", cellular = FALSE, aggregate = FALSE)
  constant <- toolHoldConstantBeyondEnd(collapseNames(input))
  getNames(constant) <- "constant"
  fadeout2050 <- constant
  fadeout2050 <- convergence(origin = constant, aim = 0, start_year = "y2020", end_year = "y2050", type = "linear")
  getNames(fadeout2050) <- "fadeout_2050"
  x <- mbind(constant, fadeout2050)

  x <- x[, findset("time"), ]

  return(list(x = x,
              weight = NULL,
              unit = "Tg Nr/yr",
              description = "Nr release through soil organic matter loss")
  )
}
