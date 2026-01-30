#' toolImportSupplyRatioAggregate
#' This function aggregates bilateral import supply ratios (amount exported from
#' country X divided by domestic supply of country Y) from country to regional
#' level, using two different weights for bilateral spatial dimensions
#' 1.1 and 1.2, where dimension 1.1 is the exporters (no weight, simple summation),
#' while dimension 1.2 are the importers, weighted by domestic supply.
#' THIS ORDER IS IMPORTANT.
#' @return return: returns region aggregated import supply ratios
#' @param x magclass object that should be aggregated
#' @param rel relation matrix containing a region mapping.
#' @param weight weight to use for dimension 1.2
#' @author David M Chen
#' @export
#'
toolImportSupplyRatioAggregate <- function(x,
                                           rel,
                                           weight) {

  rel$country <- sub("\\..*$", "", as.character(rel$country))
  rel$region  <- sub("\\..*$", "", as.character(rel$region))
  # reduce expanded bilateral mapping to unique unilateral country->region pairs
  # check for inconsistent mappings first
  counts <- vapply(split(rel$region, rel$country), function(x) length(unique(x)), integer(1))
  if (any(counts > 1)) {
    stop("rel contains countries mapping to multiple regions; inspect rel before aggregating")
  }
  # keep one row per country (unique country -> region)
  rel <- unique(rel[, c("country", "region")])


  # export is a simple summation, no weights needed
  x <- toolAggregate(x, dim = 1.1, rel = rel, from = "country", to = "region")
  # import weights are based on mass balance domestic supply
  out <- toolAggregate(x, dim = 1.2, from = "country", to = "region",
                       rel = rel, weight = weight[, getYears(x), getItems(x, dim = 3)] + 1e-10)
  return(out)
}
