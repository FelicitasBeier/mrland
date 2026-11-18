#' @title readPyrolysisCosts
#' @description Reads in raw cost parametrizations for biochar pyrolysis plant
#' set-ups from source data files.
#'
#' @param subtype Available subtypes: "REM_BCSys_Costs"
#' @return Magpie object with the dimensions "properties" and "plant_setup".
#'
#' @author Isabelle Weindl
#' @seealso [madrat::readSource()]
#' @examples
#' \dontrun{
#' x <- readSource("PyrolysisCosts", subtype = "REM_BCSys_Costs")
#' }

readPyrolysisCosts <- function(subtype = "REM_BCSys_Costs") {

  files <- c(
    REM_BCSys_Costs = "REM_BCSys_Costs.csv"
  )

  file <- toolSubtypeSelect(subtype, files)

  output <- read.magpie(file)

  getSets(output) <- c("region", "year", "properties", "plant_setup")

  return(output)
}
