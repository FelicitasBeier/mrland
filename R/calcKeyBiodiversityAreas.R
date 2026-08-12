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
#' \code{\link{readKeyBiodiversityAreas}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcKeyBiodiversityAreas", aggregate = FALSE)
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
    .alignLuDims <- function(lu, aData) {
      getYears(lu) <- NULL
      getCells(lu) <- getCells(aData)
      return(lu)
    }

    baseLandUse <- setYears(calcOutput("LanduseInitialisation", nclasses = "five",
                                       cellular = TRUE, input_magpie = TRUE,
                                       aggregate = FALSE, years = 2020),
                            "y2020")
    landUse2020 <- .alignLuDims(baseLandUse[, 2020, ], kba)

    # calculate total land area
    landArea <- dimSums(landUse2020, dim = 3)

    # urban land from LanduseInitialisation
    landUse9 <- setYears(calcOutput("LanduseInitialisation",
                                    nclasses = "nine", aggregate = FALSE,
                                    cellular = TRUE, input_magpie = TRUE,
                                    years = 2020),
                         "y2020")[, "y2020", ]
    landUse9 <- .alignLuDims(landUse9, kba)
    urbanLand2020 <- landUse9[, , "urban"]

    # make sure that kba land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand2020
    getYears(landNoUrban) <- getYears(kba)
    # compute mismatch factor
    kbaTotalLand <- dimSums(kba[, , "KBA"], dim = 3.2)
    landMismatch <- setNames(landNoUrban, NULL) / kbaTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
    # correct kba data
    kba[, , "KBA"] <- kba[, , "KBA"] * landMismatch[, , "KBA"]

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation and ESA CCI:
    kba <- toolCorrectOpenEcosystemMismatch(kba, landUse2020)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LanduseInitialisation data
      totalForest <- dimSums(landUse9[, , c("primforest", "secdforest")], dim = 3) # nolint
      primforestShr <- landUse9[, , "primforest"] / setNames(totalForest + 1e-10, NULL)
      secdforestShr <- landUse9[, , "secdforest"] / setNames(totalForest + 1e-10, NULL)
      # where LanduseInitialisation does not report forest, but we find forest land in
      # KBA data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in KBA data set
      primforest <- setNames(primforestShr, NULL) * kba[, , paste(getItems(kba, dim = 3.1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * kba[, , paste(getItems(kba, dim = 3.1), "forest", sep = ".")]

      out <- mbind(
        kba[, , c("crop", "past")],
        new.magpie(getCells(kba), getYears(kba), paste(getItems(kba, dim = 3.1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getItems(kba, dim = 3.1), "primforest", sep = ".")),
        setNames(secdforest, paste(getItems(kba, dim = 3.1), "secdforest", sep = ".")),
        new.magpie(getCells(kba), getYears(kba), paste(getItems(kba, dim = 3.1), "urban", sep = "."), fill = 0),
        kba[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      past <- new.magpie(
        cells_and_regions = getCells(kba),
        years = getYears(kba),
        names = paste(getItems(kba, dim = 3.1), "past", sep = "."),
        fill = 0
      )
      range <- kba[, , paste(getItems(kba, dim = 3.1), "past", sep = ".")]

      # separate other land into primary and secondary
      totalOther <- dimSums(landUse9[, , c("primother", "secdother")], dim = 3) # nolint
      primotherShr <- landUse9[, , "primother"] / setNames(totalOther + 1e-10, NULL)
      secdotherShr <- landUse9[, , "secdother"] / setNames(totalOther + 1e-10, NULL)
      # where LanduseInitialisation does not report other land, but we find other land in
      # KBA data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * kba[, , paste(getItems(kba, dim = 3.1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * kba[, , paste(getItems(kba, dim = 3.1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        past,
        setNames(range, paste(getItems(kba, dim = 3.1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getItems(kba, dim = 3.1), "primother", sep = ".")),
        setNames(secdother, paste(getItems(kba, dim = 3.1), "secdother", sep = "."))
      )
    }
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
