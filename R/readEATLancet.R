#' @title readEATLancet
#' @description Read in data from the EAT-Lancet Commission
#'
#' Read in data from:
#' Food in the Anthropocene: the EAT-Lancet Commission on healthy diets from sustainable food systems, Lancet 2019
#' https://doi.org/10.1016/S0140-6736(18)31788-4
#'
#'
#' @param subtype Type of EAT-Lancet data that should be read. Available types are:
#' \itemize{
#' \item \code{cons_data}: Consumption analysis ("EAT_Lancet_cons_data.csv")
#' \item \code{recommend}: Food recommendations ("EAT_Lancet_recommendations.csv")
#' }
#' @return magpie object containing EAT-Lancet Comission data
#' @importFrom data.table fread
#' @importFrom magclass complete_magpie
#' @author Isabelle Weindl, Jan Philipp Dietrich, Felicitas Beier
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' a <- readSource(type = "EATLancet", subtype = "cons_data")
#' }
#'
readEATLancet <- function(subtype) {

  if (subtype == "cons_data") {

    data <- fread("EAT_Lancet_cons_data.csv")
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

  } else if (subtype == "recommend") {

    data  <- read.csv("EAT_Lancet_recommendations_2p0.csv", sep = ",", header = TRUE,
                      stringsAsFactors = FALSE)
    mdata <- as.magpie(data, spatial = 0, temporal = 0, datacol = 2)
    getSets(mdata, fulldim = FALSE)[3] <- "target.unit.type"

  } else {
    stop("Not a valid subtype!")
  }

  return(mdata)
}
