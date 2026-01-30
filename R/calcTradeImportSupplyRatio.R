#' @title calcTradeImportSupplyRatio
#'
#' @description Calculates regional imports to supply ratios
#' often termed "Import Dependency Ratio"
#' @param magYears whether to output in magpie 5year timesteps
#' @return Self import to supply ratio
#' @author David M Chen
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTradeImportSupplyRatio()
#' }
#'
calcTradeImportSupplyRatio <- function(magYears = FALSE) {

  "!# @monitor mrland:::toolImportSupplyRatioAggregate"

  mb <- calcOutput("FAOmassbalance", aggregate = FALSE, yearly = TRUE)[, , "dm", drop = TRUE]
  tm <- calcOutput("TradeBilateralFAOHarmonized", aggregate = FALSE, yearly = TRUE)

  cyears <- intersect(getYears(tm), getYears(mb))
  citems <- intersect(getItems(tm, dim = 3), getItems(mb, dim = 3.1))

  # rename exporters so that we divide by the importers and not exporters
  getItems(tm, dim = 1.2) <- paste0(getItems(tm, dim = 1.2), "2")

  ratio <- round(tm[, cyears, citems], 6) /
    round(mb[, cyears, citems][, , "domestic_supply"], 6)
  ratio[is.na(ratio)] <- 0
  ratio[is.infinite(ratio)] <- 0
  getItems(ratio, dim = 1.2) <- gsub("[0-9]+", "", getItems(ratio, dim = 1.2))

  ratio <- collapseNames(ratio)
  ratio2 <- dimOrder(ratio, dim = 1, perm = c(2, 1))


  weight <- collapseNames(mb[, , "domestic_supply"][, cyears, citems])

  if (magYears == TRUE) {
    ratio2 <- time_interpolate(ratio2, interpolated_year = seq(1965, 1990, 5),
                               integrate_interpolated_years = TRUE,
                               extrapolation_type = "constant")
    ratio2 <- time_interpolate(ratio2, interpolated_year = seq(2025, 2150, 5),
                               integrate_interpolated_years = TRUE,
                               extrapolation_type = "constant")
    weight <- time_interpolate(weight, interpolated_year = seq(1965, 1990, 5),
                               integrate_interpolated_years = TRUE,
                               extrapolation_type = "constant")
    weight <- time_interpolate(weight, interpolated_year = seq(2025, 2150, 5),
                               integrate_interpolated_years = TRUE,
                               extrapolation_type = "constant")
    t <- magpiesets::findset("t_all")
    ratio2 <- ratio2[, t, ]
    weight <- weight [, t, ]
  }


  return(list(x = ratio2,
              weight = weight,
              unit = "ratio",
              description = "countries' import supply ratio. Imports/Domestic supply",
              aggregationFunction = toolImportSupplyRatioAggregate))
}
