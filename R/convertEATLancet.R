#' @title convertEATLancet
#' @description Convert data from the EAT Lancet Commission to be used in MAgPIE
#'
#' @param x MAgPIE object containing EAT Lancet data at mixed country-region
#' resolution
#' @param subtype Type of EAT Lancet data that should be read. Available types are:
#' \itemize{
#' \item \code{cons_data}: Consumption analysis ("EAT_Lancet_cons_data.csv")
#' \item \code{recommend}: Food recommendations ("EAT_Lancet_recommendations.csv")
#' }
#' @return EAT Lancet data as MAgPIE object at ISO country level
#' @author Isabelle Weindl, Felicitas Beier
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{
#' a <- readSource(type = "EATLancet", subtype = "cons_data")
#' }
#' @importFrom madrat getISOlist

convertEATLancet <- function(x, subtype) {

  if (subtype == "cons_data") {

    y <- toolAggregate(x, "regionmappingEATLancet.csv", weight = NULL)

  } else if (subtype == "recommend") {
    # replace NA's for min and max values
    x[, , "min"] <- ifelse(is.na(x[, , "min"]), 0, x[, , "min"])
    x[, , "max"] <- ifelse(is.na(x[, , "max"]), 999999, x[, , "max"])

    # expand to all iso countries
    iso     <- getISOlist(type = "all")
    mapping <- cbind(rep("GLO", length(iso)), as.character(iso))
    y       <- toolAggregate(x, rel = mapping, weight = NULL)

  } else {

    stop("Not a valid subtype! Please select cons_data or recommend in readEATLancet function.")

  }

  return(y)
}
