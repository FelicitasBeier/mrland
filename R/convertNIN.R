#' Convert data from the NIN Lancet Comission
#'
#' Convert data from the NIN Lancet Comission to ISO country level.
#'
#'
#' @param x MAgPIE object containing NIN Lancet data at mixed country-region
#' resolution
#' @param subtype Type of NIN Lancet data that should be read. Available types are:
#' \itemize{
#' \item \code{cons_data}: Consumption analysis ("NIN_Lancet_cons_data.csv")
#' \item \code{recommend}: Food recommendations ("NIN_recommendations.csv")
#' }
#' @return NIN Lancet data as MAgPIE object at ISO country level
#' @author Isabelle Weindl
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#'
#' \dontrun{
#' a <- readSource(type="NIN",subtype="cons_data")
#' }
#' @importFrom madrat getISOlist

convertNIN <- function(x,subtype) {

if(subtype == "cons_data"){

    y <- toolAggregate(x, "regionmappingEATLancet.csv", weight=NULL )

} else if (subtype == "recommend") {

  iso <- getISOlist(type = "all")
  mapping <- cbind(rep("GLO",length(iso)),as.character(iso))
  y <- toolAggregate(x, rel = mapping, weight=NULL )

} else {

  stop("Not a valid subtype!")

}

  return(y)
}









