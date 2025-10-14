#' @title calcTradeSelfSuff
#'
#' @description Calculates regional self sufficiences from FAO data as
#' production/domestic_supply.
#'
#' @return Self sufficiences
#' @author Ulrich Kreidenweis
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link[mrcommons]{calcFAOmassbalance}}
#' @examples
#' \dontrun{
#' a <- calcTradeSelfSuff()
#' }
#'
calcTradeSelfSuff <- function() {

  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)
  # add missing products - only scp at the moment
  massbalance[, , "scp"] <- 0
  massbalance <- massbalance[, , findset("k_trade")]

  selfSuff <- massbalance[, , "production.dm"] / massbalance[, , "domestic_supply.dm"]
  selfSuff <- collapseNames(selfSuff)
  selfSuff[is.nan(selfSuff)] <- 0
  selfSuff[selfSuff == Inf] <- 1
  selfSuff[, , "scp"] <- 1

  weight <- massbalance[, , "domestic_supply.dm"]
  weight <- collapseNames(weight)
  weight[is.nan(weight)] <- 0
  weight[, , "scp"] <- 1

  out <- toolHoldConstantBeyondEnd(selfSuff)
  weight <- toolHoldConstantBeyondEnd(weight)
  weight <- weight + 10^-10

  return(list(x = out,
              weight = weight,
              unit = "ratio",
              description = "countries' self sufficiencies in agricultural production. Production/Domestic supply"))
}
