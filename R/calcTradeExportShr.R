#' @title calcTradeExportShr
#' @description Calculate  export shares of countries compared to total export. 
#'              This is based on export values from FAOSTAT. 
#'              Function calculates this based on average
#'              values of the specified years.
#'
#' @return export shares as MAgPIE object
#' @author Ulrich Kreidenweis, Xiaoxi Wang
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTradeExportShr()
#' }
#'
calcTradeExportShr <- function() {

  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)
  massbalance <- massbalance[, , findset("k_trade")]

  # export share is the regions' share in total exports. 
  # So the sum per commodity over all regions is 1
  netexp <- massbalance[, , "export"][, , "dm"] - massbalance[, , "import"][, , "dm"]
  netexp[netexp < 0] <- 0
  netexp[, , "scp"] <- 1

  expGlo <- dimSums(netexp, dim = 1)
  expShr <- netexp / expGlo
  expShr <- collapseNames(expShr)
  expShr[is.nan(expShr)] <- 0

  expShr <- toolHoldConstantBeyondEnd(expShr)

  return(list(x = expShr,
              weight = NULL,
              unit = "share",
              description = "share of export of individual counties in total global exports"))
}
