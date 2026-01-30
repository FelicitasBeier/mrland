#' @title calcTradeBilateralBalanceFlow
#'
#' @description Calculates balanceflows for trade and total mass balance
#' to match historical production
#' @param balanceflow "trade" or "total"
#' @return Self import to dupply ratio
#' @author David M Chen
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTrademportSupplyStdDev()
#' }
#' @importFrom dplyr %>% inner_join
#' @importFrom rlang .data :=

calcTradeBilateralBalanceFlow <- function(balanceflow = "trade") {

  tm <- calcOutput("TradeBilateralFAOHarmonized", aggregate = FALSE, yearly = TRUE)
  mb <- calcOutput("FAOmassbalance", aggregate = FALSE)[, , "dm", drop = TRUE]

  # tm and mb imports should match
  cyears <- intersect(getYears(tm), getYears(mb))
  citems <- intersect(getItems(mb, dim = 3.1), getItems(tm, dim = 3))

  if (any(round(dimSums(tm[, cyears, citems], dim = 1.2) -
                  mb[, cyears, citems][, , "import"], 3) != 0)) {
    warning("FAO trade matrix imports and mass balance imports don't match, \n
           something wrong in the TradeBilateralFAOHarmonized fucntion")
  }

  tmx <- dimSums(tm, dim = 1.1)
  mbx <- mb[, , "export", drop = TRUE]
  out <- mbx[, cyears, citems] - tmx[, cyears, citems]
  xbf <- out

  if (balanceflow == "total") {
    # Exports and imports also mismatch in the mass balance so we need to then
    # account for full mass balance prod - domestic supply + exports - imports
    # balance flow

    out <- collapseNames((mb[, cyears, citems][, , "production"] -
                            mb[, cyears, citems][, , "domestic_supply"] -
                            xbf[, cyears, citems] -
                            tmx[, cyears, citems] +
                            mb[, cyears, citems][, , "import"]))
  }

  # make years after history 0
  out <- toolHoldConstantBeyondEnd(out)
  # fading out the balanceflow until 2030.
  out <- convergence(origin = out, aim = 0,
                     start_year = "y2020", end_year = "y2030",
                     type = "s")

  return(list(x = out,
              weight = NULL,
              unit = "Mt",
              description = "Balanceflow for bilateral trade matrix imports not matching exports"))

}
