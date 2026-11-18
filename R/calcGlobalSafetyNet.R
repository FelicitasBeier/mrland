#' @title calcGlobalSafetyNet
#'
#' @description Returns unprotected land area (Mha) within the Global Safety Net (Dinerstein et al. 2020).

#' @param maginput Whether data should be transformed (based on LUH3 data) to match land use types used in MAgPIE.
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
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{readDinerstein2020}}, \code{\link{toolMagpieLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("GlobalSafetyNet", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcGlobalSafetyNet <- function(maginput = TRUE, nclasses = "seven", cells = "lpjcell") {
  gsn <- mbind(
    readSource("Dinerstein2020", subtype = "GSN:distinct_species_assemblages", convert = "onlycorrect"),
    readSource("Dinerstein2020", subtype = "GSN:rare_phenomena", convert = "onlycorrect"),
    readSource("Dinerstein2020", subtype = "GSN:areas_of_intactness", convert = "onlycorrect"),
    readSource("Dinerstein2020", subtype = "GSN:climate_stabilisation_tier1", convert = "onlycorrect"),
    readSource("Dinerstein2020", subtype = "GSN:climate_stabilisation_tier2", convert = "onlycorrect")
  )

  if (maginput == TRUE) {
    out <- toolMagpieLandClasses(gsn, nclasses)
  } else {
    out <- gsn
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = paste(
      "Unprotected land area of the Global Safety Net (Dinerstein et al. 2020)."
    ),
    isocountries = FALSE
  ))
}
