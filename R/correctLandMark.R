#' @title correctLandMark
#' @description correct data of the LandMark IPLC lands.
#' @return magpie object on cellular level
#' @param x magpie object provided by the read function
#' @author Patrick v. Jeetze
#' @seealso
#'   \code{\link{readLandMark}}
#' @examples
#' \dontrun{
#' readSource("LandMark", convert = "onlycorrect")
#' }
#'
#' @importFrom madrat toolConditionalReplace

correctLandMark <- function(x) {
  x <- toolConditionalReplace(x, conditions = "is.na()", replaceby = 0)

  return(x)
}
