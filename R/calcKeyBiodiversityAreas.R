#' @title calcKeyBiodiversityAreas
#'
#' @description Returns unprotected land area (Mha) within Key Biodiversity Areas.

#' @param maginput Whether data should be transformed (based on luh3 data) to match land use types used in MAgPIE.
#' @param nclasses If \code{magpie_input = TRUE}. Options are either "seven" or "nine". Note that by default,
#' the protected area is reported for urban land and forestry is zero.
#' \itemize{
#' \item "seven" separates primary and secondary forest and includes
#' "crop", "past", "forestry", "primforest", "secdforest", "urban" and "other"
#' \item "nine" adds the separation of pasture and rangelands, as well as a
#' differentiation of primary and secondary non-forest vegetation and therefore returns
#' "crop", "past", "range", "forestry", "primforest", "secdforest", "urban", "primother" and "secdother"
#' }
#' @param unprotected if TRUE only KBA land that is currently unprotected is returned
#' @param cells magpiecell (59199 cells) or lpjcell (67420 cells)
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{readKeyBiodiversityAreas}}, \code{\link{toolMagpieLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("KeyBiodiversityAreas", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcKeyBiodiversityAreas <- function(maginput = TRUE, unprotected = TRUE,
                                     nclasses = "seven", cells = "lpjcell") {
  if (unprotected) {
    kba <- readSource("KeyBiodiversityAreas", subtype = "unprotected", convert = "onlycorrect")
  } else {
    kba <- readSource("KeyBiodiversityAreas", subtype = "all", convert = "onlycorrect")
  }

  if (maginput) {
    out <- toolMagpieLandClasses(kba, nclasses)
  } else {
    out <- kba
  }

  if (cells == "magpiecell") {
    out <- toolCoord2Isocell(out)
  } else if (cells != "lpjcell") {
    stop("Please specify cells argument")
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = paste(
      "Unprotected land area in Key Biodiversity Areas"
    ),
    isocountries = FALSE
  ))
}
