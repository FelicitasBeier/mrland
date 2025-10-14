#' @title Calculate global under-/overproduction
#'
#' @description Calculate the difference between the global production and
#'              global domestic_supply (corrective balance flow).
#'              The difference is the result of imports not equaling exports,
#'              and because storage is not considered. The calculated
#'              DomesticBalanceflow assures that production matches domestic_supply. The
#'              goods come from nowhere and go to nowhere. The numbers are usually decreased
#'              linearly and become zero in 2050.
#'
#' @return global domestic balanceflows as MAgPIE object
#' @author Ulrich Kreidenweis, Xiaoxi Wang
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTradeBalanceflow()
#' }
#'
calcTradeBalanceflow <- function() {

  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)
  kTrade      <- findset("k_trade")
  massbalance <- massbalance[, , kTrade]

  x <- dimSums(massbalance[, , c("production.dm", "domestic_supply.dm")], dim = 1)
  balanceflow <- collapseNames(x[, , "production.dm"]) - collapseNames(x[, , "domestic_supply.dm"])

  out <- toolHoldConstantBeyondEnd(balanceflow)
  # fading out the balanceflow until 2030.
  out <- convergence(origin = out, aim = 0,
                     start_year = "y2020", end_year = "y2030",
                     type = "s")

  return(list(x = out,
              weight = NULL,
              unit = "mio. ton dm",
              description = "Balanceflow to match global production and domestic supply"))
}
