#' Convert FRA2015Doc data
#'
#'
#' @param x MAgPIE object containing original values coming from read function
#' @param subtype The data table type, e.g.: forest_area
#' @return Data as MAgPIE object
#' @author Abhijeet Mishra
#' @seealso \code{\link{readFRA2015Doc}}, \code{\link[madrat]{readSource}},
#' @examples
#' \dontrun{
#' a <- readSource("FRA2015Doc", "forest_area", convert = TRUE)
#' }
#' @importFrom magclass getItems
#' @importFrom madrat toolCountryFill
#'


convertFRA2015Doc <- function(x, subtype) {

  unreported <- FALSE
  missing    <- FALSE

 if (!is.null(x)) {
   fill <- mean(x, na.rm = TRUE)
   x[is.na(x)] <- 0

   ## missing data
   for (i in getItems(x, dim = 1.1)) {
     if (all(x[i, , ] == 0)) {
       x[i, , ] <- fill
       unreported <- TRUE
       }
     if (any(x[i, , ] == 0)) {
       x[i, , ] <- max(x[i, , ])
       missing <- TRUE
       }
   }
   if (unreported) message("Countries with no data in FRA 2015 have been given the mean value.\n")
   if (missing) message("Countries with missing data in FRA 2015 have been given the maximum reported value.")

   x <- toolCountryFill(x, fill = 0, verbosity = 2)
  } else {
return(x)
}
}
