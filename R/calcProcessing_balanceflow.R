#' @title calcProcessing_balanceflow
#'
#' @description Calculates the difference between
#' what is recorded in FAOSTAT production of secondary
#' products, and what is calculated from our global
#' processing conversion factors in order to meet FAOSTAT
#' production. Fades out by 2050.
#'
#' @return magpie object of secondary product balanceflow
#' @author Benjamin Bodirsky, David M Chen
#' @seealso \code{\link[mrland]{calcProcessing_conversion_factors}}

#' @examples
#' \dontrun{
#' a <- calcOutput("Processing_balanceflow",
#'                 aggregate = FALSE)
#' }

calcProcessing_balanceflow <- function() { # nolint: object_name_linter.
  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)

  # add cotton fibres, cause its not in the massabalance calculations
  conv <- calcOutput("Processing_conversion_factors", aggregate = FALSE)[, , "ginning"]
  massbalance[, , "fibres.production_estimated.dm"] <-
    collapseNames(conv[, getYears(massbalance), "ginning.fibres.cottn_pro"] *
                    massbalance[, , "cottn_pro.production.dm"])

  ksd <- findset("ksd")
  ksdNoscp <- setdiff(ksd, "scp")
  balanceflow <- dimSums(massbalance[, , "production"][, , "dm"][, , ksdNoscp], dim = c(3.2, 3.3)) -
    dimSums(massbalance[, , "production_estimated"][, , "dm"][, , ksdNoscp], dim = c(3.2, 3.3))
  balanceflow <- add_columns(balanceflow, addnm = "scp", dim = 3.1)
  balanceflow[, , "scp"] <- 0

  out <- toolHoldConstantBeyondEnd(balanceflow)
  # fading out the balanceflow until 2050.
  out <- convergence(origin = out, aim = 0, start_year = "y2015", end_year = "y2050", type = "s")

  return(list(x = out, weight = NULL, unit = "t DM",
              description = paste0("Balanceflow to balance out the assumption of uniform ",
                                   "conversion factors worldwide but match FAO.")))
}
