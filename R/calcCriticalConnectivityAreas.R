#' @title calcCriticalConnectivityAreas
#'
#' @description Returns unprotected land area (Mha) within Critical Connectivit Areas
#' as given in Brennan et al. (2022).

#' @param maginput Whether data should be transformed (based on LanduseInitialisation data) to match
#' land use types used in MAgPIE.
#' @param nclasses If \code{magpie_input = TRUE}. Options are either "seven" or "nine". Note that by default,
#' the protected area is reported for urban land and forestry is zero.
#' \itemize{
#' \item "seven" separates primary and secondary forest and includes
#' "crop", "past", "forestry", "primforest", "secdforest", "urban" and "other"
#' \item "nine" adds the separation of pasture and rangelands, as well as a
#' differentiation of primary and secondary non-forest vegetation and therefore returns
#' "crop", "past", "range", "forestry", "primforest", "secdforest", "urban", "primother" and "secdother"
#' }
#' @param cells (deprecated) always lpjcell (67420 cells)
#' @param mask Whether Key Biodiversity Areas ("KBA") or Global Safety Net and Key Biodiversity Areas
#' ("KBA_GSN") are masked. This switch is useful for complementary scenario building.
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{readBrennan2022}}, \code{\link{toolMagpieLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcCriticalConnectivityAreas", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcCriticalConnectivityAreas <- function(maginput = TRUE, nclasses = "seven",
                                          cells = "lpjcell", mask = "KBA_GSN") {
  if (mask == "KBA") {
    cca <- readSource("Brennan2022", subtype = "KBA_masked", convert = "onlycorrect")
  } else if (mask == "KBA_GSN") {
    cca <- readSource("Brennan2022", subtype = "KBA_GSN_masked", convert = "onlycorrect")
  } else {
    stop("Option specified for argument 'mask' does not exist.")
  }

  if (maginput == TRUE) {
    out <- toolMagpieLandClasses(cca, nclasses)
  } else {
    out <- cca
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = paste(
      "Unprotected land area in Critical Connectivity Areas (Brennan et al. 2022)."
    ),
    isocountries = FALSE
  ))
}
