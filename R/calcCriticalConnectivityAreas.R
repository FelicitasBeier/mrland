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
#' \code{\link{readBrennan2022}}
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
    .alignLuDims <- function(lu, aData) {
      getYears(lu) <- NULL
      getCells(lu) <- getCells(aData)
      return(lu)
    }

    baseLandUse <- calcOutput("LanduseInitialisation", nclasses = "five", cellular = TRUE,
                              input_magpie = TRUE, aggregate = FALSE, years = 2020)
    landUse2020 <- .alignLuDims(baseLandUse[, 2020, ], cca)

    # calculate total land area
    landArea <- dimSums(landUse2020, dim = 3)

    # urban land from LanduseInitialisation
    landUse9 <- calcOutput("LanduseInitialisation",
      nclasses = "nine", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE,
      years = 2020
    )[, "y2020", ]
    landUse9 <- .alignLuDims(landUse9, cca)
    urbanLand2020 <- landUse9[, , "urban"]

    # make sure that cca land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand2020
    getYears(landNoUrban) <- getYears(cca)
    # compute mismatch factor
    ccaTotalLand <- dimSums(cca[, , "CCA"], dim = 3.2)
    landMismatch <- setNames(landNoUrban, NULL) / ccaTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
    # correct cca data
    cca[, , "CCA"] <- cca[, , "CCA"] * landMismatch[, , "CCA"]

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation and ESA CCI:
    cca <- toolCorrectOpenEcosystemMismatch(cca, landUse2020)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LanduseInitialisation data
      totalForest <- dimSums(landUse9[, , c("primforest", "secdforest")], dim = 3) # nolint
      primforestShr <- landUse9[, , "primforest"] / setNames(totalForest + 1e-10, NULL)
      secdforestShr <- landUse9[, , "secdforest"] / setNames(totalForest + 1e-10, NULL)
      # where LanduseInitialisation does not report forest, but we find forest land in
      # CCA data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in CCA data set
      primforest <- setNames(primforestShr, NULL) * cca[, , paste(getItems(cca, dim = 3.1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * cca[, , paste(getItems(cca, dim = 3.1), "forest", sep = ".")]

      out <- mbind(
        cca[, , c("crop", "past")],
        new.magpie(getCells(cca), getYears(cca), paste(getItems(cca, dim = 3.1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getItems(cca, dim = 3.1), "primforest", sep = ".")),
        setNames(secdforest, paste(getItems(cca, dim = 3.1), "secdforest", sep = ".")),
        new.magpie(getCells(cca), getYears(cca), paste(getItems(cca, dim = 3.1), "urban", sep = "."), fill = 0),
        cca[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      past <- new.magpie(
        cells_and_regions = getCells(cca),
        years = getYears(cca),
        names = paste(getItems(cca, dim = 3.1), "past", sep = "."),
        fill = 0
      )
      range <- cca[, , paste(getItems(cca, dim = 3.1), "past", sep = ".")]

      # separate other land into primary and secondary
      totalOther <- dimSums(landUse9[, , c("primother", "secdother")], dim = 3) # nolint
      primotherShr <- landUse9[, , "primother"] / setNames(totalOther + 1e-10, NULL)
      secdotherShr <- landUse9[, , "secdother"] / setNames(totalOther + 1e-10, NULL)
      # where LanduseInitialisation does not report other land, but we find other land in
      # CCA data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * cca[, , paste(getItems(cca, dim = 3.1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * cca[, , paste(getItems(cca, dim = 3.1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        past,
        setNames(range, paste(getItems(cca, dim = 3.1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getItems(cca, dim = 3.1), "primother", sep = ".")),
        setNames(secdother, paste(getItems(cca, dim = 3.1), "secdother", sep = "."))
      )
    }
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
