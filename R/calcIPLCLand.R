#' @title calcIPLCLand
#'
#' @description Returns unprotected land area (Mha) covered by Indigenous Peoples' and Local Community Lands.

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
#' @param datasource Currently only \code{"LandMark"}
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{readLandMark}}, \code{\link{toolMagpieLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcIPLCLand", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcIPLCLand <- function(maginput = TRUE, nclasses = "seven", datasource = "LandMark") {
  if (datasource == "LandMark") {
    iplc <- mbind(
      readSource("LandMark", subtype = "delineated", convert = "onlycorrect"),
      readSource("LandMark", subtype = "indicative", convert = "onlycorrect")
    )
    iplcAll <- dimSums(iplc, dim = 3.1)
    iplcAll <- add_dimension(iplcAll, dim = 3.1, nm = "LandMark_IPLC_all")
    iplc <- mbind(iplc, iplcAll)
  } else {
    stop("Please select a valid data source")
  }


  if (maginput == TRUE) {
    out <- toolMagpieLandClasses(iplc, nclasses)
  } else {
    out <- iplc
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = paste(
      "Unprotected land area covered by Indigenous Peoples' and Local Community Lands",
      paste0("(source: '", datasource, "')")
    ),
    isocountries = FALSE
  ))
}
