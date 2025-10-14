#' Read in data from the NIN recommendations
#'
#' @param subtype Type of NIN data that should be read. Available types are:
#' \itemize{
#' \item \code{cons_data}: Consumption analysis ("NIN_cons_data.csv")
#' }
#' @return magpie object containing NIN data
#' @author Isabelle Weindl, Jan Philipp Dietrich
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' a <- readSource(type = "NIN", subtype = "cons_data")
#' }
#'
readNIN <- function(subtype) {
  if (subtype == "cons_data") {
    data <- data.table::fread("NIN_cons_data.csv")
    data$SSP_scn <- NULL # remove irrelevant column
    data <- data[!(data$region %in% c("all-r", "HIC", "UMC", "LMC", "LIC")), ] # remove aggregated regions

    # only keep the absolute values (abs) and remove the relative measures
    # (absolute changes relative to 2010 (chg_2010), percentage changes relative to 2010 (pct_2010))
    data <- data[data$measure == "abs", ]
    data$measure <- NULL

    # re-structure the dataframe such that for all scenarios and units the respective data sets have the same size
    data <- data[data$food_group != "total", ]

    # Fill gaps with 0
    mdata <- as.magpie(data, tidy = TRUE)
    mdata <- complete_magpie(mdata)
    mdata[is.na(mdata)] <- 0
  } else {
    stop("Not a valid subtype!")
  }

  return(mdata)
}
