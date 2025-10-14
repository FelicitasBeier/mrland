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
    luh3 <- calcOutput("LUH3",
      landuseTypes = "LUH3", aggregate = FALSE,
      cellular = TRUE,
      yrs = 2015
    )
    getYears(luh3) <- NULL
    getCells(luh3) <- getCells(kba)

    # calculate total land area
    landArea <- dimSums(luh3, dim = 3)

    # urban land
    urbanLand <- calcOutput("UrbanLandFuture",
      subtype = "LUH3", aggregate = FALSE,
      timestep = "5year", cells = "lpjcell"
    )

    # make sure that kba land is not greater than total land area minus urban area
    landNoUrban <- setYears(landArea, "y2020") - setCells(urbanLand[, "y2020", "SSP2"], getCells(landArea))
    getYears(landNoUrban) <- getYears(kba)
    # compute mismatch factor
    kbaTotalLand <- dimSums(kba[, , "KBA"], dim = 3.2)
    landMismatch <- setNames(landNoUrban, NULL) / kbaTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
    # correct kba data
    kba[, , "KBA"] <- kba[, , "KBA"] * landMismatch[, , "KBA"]

    if (nclasses %in% c("seven", "nine")) {

      # differentiate primary and secondary forest based on luh3 data
      totForestLUH <- dimSums(luh3[, , c("primf", "secdf")], dim = 3) # nolint
      primforestShr <- luh3[, , "primf"] / setNames(totForestLUH + 1e-10, NULL)
      secdforestShr <- luh3[, , "secdf"] / setNames(totForestLUH + 1e-10, NULL)
      # where luh2 does not report forest, but we find forest land in
      # KBA data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in KBA data set
      primforest <- setNames(primforestShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "forest", sep = ".")]

      out <- mbind(
        kba[, , c("crop", "past")],
        new.magpie(getCells(kba), getYears(kba), paste(getNames(kba, dim = 1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getNames(kba, dim = 1), "primforest", sep = ".")),
        setNames(secdforest, paste(getNames(kba, dim = 1), "secdforest", sep = ".")),
        new.magpie(getCells(kba), getYears(kba), paste(getNames(kba, dim = 1), "urban", sep = "."), fill = 0),
        kba[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {

      # separate pasture into pasture and rangeland
      totGrassLUH <- dimSums(luh3[, , c("pastr", "range")], dim = 3) # nolint
      pastShr <- luh3[, , "pastr"] / setNames(totGrassLUH + 1e-10, NULL)
      rangeShr <- luh3[, , "range"] / setNames(totGrassLUH + 1e-10, NULL)
      # where luh2 does not report grassland, but we find grassland in
      # KBA data, set share of rangeland to 1
      rangeShr[pastShr == 0 & rangeShr == 0] <- 1
      # multiply shares of pasture and rangeland with pasture in KBA data
      past <- setNames(pastShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "past", sep = ".")]
      range <- setNames(rangeShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "past", sep = ".")]

      # separate other land into primary and secondary
      totOtherLUH <- dimSums(luh3[, , c("primn", "secdn")], dim = 3) # nolint
      primotherShr <- luh3[, , "primn"] / setNames(totOtherLUH + 1e-10, NULL)
      secdotherShr <- luh3[, , "secdn"] / setNames(totOtherLUH + 1e-10, NULL)
      # where luh2 does not report other land, but we find other land in
      # KBA data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * kba[, , paste(getNames(kba, dim = 1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        setNames(past, paste(getNames(kba, dim = 1), "past", sep = ".")),
        setNames(range, paste(getNames(kba, dim = 1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getNames(kba, dim = 1), "primother", sep = ".")),
        setNames(secdother, paste(getNames(kba, dim = 1), "secdother", sep = "."))
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
