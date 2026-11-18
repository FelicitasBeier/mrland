#' @title calcIrrecoverableCarbonLand
#'
#' @description Returns unprotected land area (Mha) that covers 50 %, 75% and
#' 99% of irrecoverable carbon as defined in Noon et al (2022).

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
#' @param cells magpiecell (59199 cells) or lpjcell (67420 cells)
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{readNoon2022}}, \code{\link{toolMagpieLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("IrrecoverableCarbonLand", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcIrrecoverableCarbonLand <- function(maginput = TRUE, nclasses = "seven",
                                        cells = "lpjcell") {
  ic <- mbind(
    readSource("Noon2022", subtype = "land:IrrC_50pc", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_75pc", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_95pc", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_99pc", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_75pc_30by30", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_95pc_30by30", convert = "onlycorrect"),
    readSource("Noon2022", subtype = "land:IrrC_99pc_30by30", convert = "onlycorrect")
  )

  if (maginput == TRUE) {
    out <- toolMagpieLandClasses(ic, nclasses)
  } else {
    out <- ic
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
      "Unprotected land area that contains",
      "50 %, 75 %, 95 % and 99 %",
      "of irrecoverable carbon as defined in Noon et al (2022)."
    ),
    isocountries = FALSE
  ))
}
