#' @title correctS4Nproject_input
#' @description corrects IMAGE inputs of total bioenergy (1st gen, 2nd gen and residues) demand and co2 prices
#' 
#' @param x magpie object
#' 
#' @return magpie object at country-level resolution
#' @author Felicitas Beier
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' 
#' \dontrun{ a <- readSource("S4Nproject_input", aggregate=FALSE)
#' }

correctS4Nproject_input <- function(x) {
 
  # fill missing countries
  x <- toolCountryFill(x, fill=0)
  # fill missing years
  x <- toolFillYears(x,paste0("y", seq(2005, 2100,by=5)))
  
  return(x)
}
