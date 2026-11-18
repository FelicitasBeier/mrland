#' @title calcProtectedAreaBaseline
#'
#' @description Returns protected land area (Mha) in terms of cropland, pasture, forest and other
#' land between 1995 and 2020.
#'
#' @param magpie_input Whether data should be transformed (based on LUH3 data) to match land use types used in MAgPIE.
#' @param nclasses If \code{magpie_input = TRUE}. Options are either "seven" or "nine". Note that by default,
#' the protected area is reported for urban land and forestry is zero.
#' \itemize{
#' \item "seven" separates primary and secondary forest and includes
#' "crop", "past", "forestry", "primforest", "secdforest", "urban" and "other"
#' \item "nine" adds the separation of pasture and rangelands, as well as a
#' differentiation of primary and secondary non-forest vegetation and therefore returns
#' "crop", "past", "range", "forestry", "primforest", "secdforest", "urban", "primother" and "secdother"
#' }
#' @param cells (deprecated) only lpjcell (67420 cells)
#'
#' @return List with a magpie object
#' @author Patrick v. Jeetze
#' @seealso
#' \code{\link{calcProtectArea}}, \code{\link{toolSplitLandClasses}}
#'
#' @examples
#' \dontrun{
#' calcOutput("ProtectedAreaBaseline", aggregate = FALSE)
#' }
#'
#' @importFrom mstools toolCoord2Isocell
#'
calcProtectedAreaBaseline <- function(magpie_input = TRUE, nclasses = "seven", # nolint
                                      cells = "lpjcell") {
  PABaseline <- readSource("ProtectedAreaBaseline", convert = "onlycorrect") # nolint

  if (magpie_input == TRUE) {

    landUse9 <- calcOutput("LanduseInitialisation",
      nclasses = "nine", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE,
      years = c(1995, 2000, 2005, 2010, 2015, 2020)
    )

    # extend the data set to all time steps provided in the protected area data
    # i.e. use the data from the year 2015 for the year 2020.
    landUse9 <- setCells(landUse9, getCells(PABaseline))

    # calculate total land area
    landArea <- dimSums(landUse9[, "y1995", ], dim = 3)

    # urban land
    urbanLand <- calcOutput("UrbanLandFuture",
      subtype = "LUH3", aggregate = FALSE,
      timestep = "5year", cells = "lpjcell"
    )

    # make sure that protected area is not greater than total land area minus urban area
    totPABase <- dimSums(PABaseline, dim = 3.2)
    landNoUrban <- landArea - setCells(urbanLand[, getYears(PABaseline), "SSP2"], getCells(PABaseline))
    getYears(landNoUrban) <- getYears(PABaseline)
    # compute mismatch factor
    landMismatch <- setNames(landNoUrban, NULL) / totPABase
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
    # correct WDPA data
    PABaseline <- PABaseline * landMismatch # nolint

    ### Derive area of remaining protected area categories IV-V-VI
    # account for small mismatches in categories
    catMismatch <- PABaseline[, , "WDPA"] < PABaseline[, , "WDPA_I-II-III"]
    PABaseline[catMismatch] <- PABaseline[, , "WDPA"][catMismatch] # nolint
    PABaseline <- mbind( # nolint
      PABaseline,
      setNames(
        PABaseline[, , "WDPA"] - PABaseline[, , "WDPA_I-II-III"],
        sub("WDPA", "WDPA_IV-V-VI", getNames(PABaseline[, , "WDPA"]))
      )
    )

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation (LUH) and ESA CCI:
    luIni <- calcOutput("LanduseInitialisation",
      nclasses = "seven", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE
    )[, getYears(PABaseline), ]
    luIni <- setCells(luIni, getCells(PABaseline))
    PABaseline <- toolCorrectOpenEcosystemMismatch(PABaseline, luIni) # nolint

    out <- toolSplitLandClasses(PABaseline, landUse9, nclasses)
  } else {
    out <- PABaseline
  }

  return(list(
    x = out,
    weight = NULL,
    unit = "Mha",
    description = "Protected land area in terms of cropland, pasture, forest and other land between 1995 and 2020.",
    isocountries = FALSE
  ))
}
