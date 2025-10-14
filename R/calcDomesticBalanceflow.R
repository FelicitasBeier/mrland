#' @title calcDomesticBalanceflow
#' @description calculate the corrective balance flow that is needed
#' because global production doesn't equal global supply
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Benjamin Leon Bodirsky
#' @seealso
#' [calcFAOmassbalance()]
#' @examples
#' \dontrun{
#' calcOutput("DomesticBalanceflow")
#' }

#' @importFrom magclass dimSums convergence

calcDomesticBalanceflow <- function() {
  mb <- calcOutput("FAOmassbalance", aggregate = F)
  processing <- c("milling", "refining", "extracting", "fermentation", "distilling")

  mb2 <- mb[, , c("food", "feed", "seed", "waste", "other_util", "bioenergy")][, , "dm"]
  mb2 <- add_columns(x = mb2, dim = 3.2, addnm = "processed")
  mb2[, , "processed"] <- dimSums(mb[, , processing][, , "dm"], dim = 3.2)

  balanceflow <- collapseNames(mb[, , "domestic_supply"][, , "dm"]) - collapseNames(dimSums(mb2[, , ], dim = c(3.2)))

  out <- toolHoldConstantBeyondEnd(balanceflow)
  # fading out the balanceflow until 2050.
  out <- convergence(origin = out, aim = 0, start_year = "y2020", end_year = "y2050", type = "s")

  return(list(
    x = out,
    weight = NULL,
    unit = "mio. ton dm",
    description = "Balanceflow to domestic supply by different demand categories"
  ))
}
