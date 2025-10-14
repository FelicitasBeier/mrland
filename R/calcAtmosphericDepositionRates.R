#' @title calcAtmosphericDepositionRates
#' @description Conputes Atmospheric (nitrogen) deposition rates per area on different land-use types.
#'
#' @param cellular TRUE for results on 0.5 degree grid.
#' @param cells            magpiecell (59199 cells) or lpjcell (67420 cells)
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Benjamin Leon Bodirsky
#'
#' @seealso
#' \code{\link[mrcommons]{calcAtmosphericDeposition}},
#' \code{\link[mrcommons]{calcNitrogenBudgetCropland}}
#' @examples
#' \dontrun{
#' calcOutput("AtmosphericDepositionRates")
#' }
#'
#' @importFrom magclass nregions

calcAtmosphericDepositionRates <- function(cellular = FALSE, cells = "lpjcell") {

  # dep <- calcOutput("AtmosphericDeposition", datasource = "ACCMIP",
  #                   scenario = c("rcp26", "rcp45", "rcp85"), aggregate = FALSE)
  dep     <- calcOutput("AtmosphericDeposition", datasource = "Nsurplus2",
                        cellular = cellular, cells = cells, aggregate = FALSE)
  dep     <- dimSums(dep, dim = c(3.3, 3.4))
  luhdata <- calcOutput("LanduseInitialisation", cellular = cellular,
                        cells = cells, aggregate = FALSE)
  dep     <- toolHoldConstantBeyondEnd(dep)
  weight  <- toolHoldConstantBeyondEnd(luhdata)
  dep     <- dep * (weight > 10^-7)

  out <- dep / weight
  out[is.na(out)] <- 0
  out[is.infinite(out)] <- 0

  weight <- weight + 10^(-10)

  return(list(x = out,
              weight = weight,
              unit = "Mt Nr / Mha",
              min = 0,
              max = 2000,
              isocountries = (!cellular & (nregions(out) != 1)),
              description = "Atmospheric deposition per ha on different land types."))
}
