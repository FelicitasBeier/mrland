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
#' \code{\link{readDinerstein2020}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcGlobalSafetyNet", aggregate = FALSE)
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
    .alignLuDims <- function(lu, aData) {
      getYears(lu) <- NULL
      getCells(lu) <- getCells(aData)
      return(lu)
    }

    baseLandUse <- calcOutput("LanduseInitialisation", nclasses = "five", cellular = TRUE,
                              input_magpie = TRUE, aggregate = FALSE, years = 2020)
    landUse2020 <- .alignLuDims(baseLandUse[, 2020, ], gsn)

    # calculate total land area
    landArea <- dimSums(landUse2020, dim = 3)

    # urban land from LanduseInitialisation
    landUse9 <- calcOutput("LanduseInitialisation",
      nclasses = "nine", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE,
      years = 2020
    )[, "y2020", ]
    landUse9 <- .alignLuDims(landUse9, gsn)
    urbanLand2020 <- landUse9[, , "urban"]

    # make sure that GSN land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand2020
    getYears(landNoUrban) <- getYears(gsn)

    # compute mismatch factor
    gsnTotalLand <- mbind(
      dimSums(gsn[, , "GSN_DSA"], dim = 3.2),
      dimSums(gsn[, , "GSN_RarePhen"], dim = 3.2),
      dimSums(gsn[, , "GSN_AreaIntct"], dim = 3.2),
      dimSums(gsn[, , "GSN_ClimTier1"], dim = 3.2),
      dimSums(gsn[, , "GSN_ClimTier2"], dim = 3.2)
    )
    landMismatch <- setNames(landNoUrban, NULL) / gsnTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)

    # correct GSN data
    gsn[, , "GSN_DSA"] <- gsn[, , "GSN_DSA"] * landMismatch[, , "GSN_DSA"]
    gsn[, , "GSN_RarePhen"] <- gsn[, , "GSN_RarePhen"] * landMismatch[, , "GSN_RarePhen"]
    gsn[, , "GSN_AreaIntct"] <- gsn[, , "GSN_AreaIntct"] * landMismatch[, , "GSN_AreaIntct"]
    gsn[, , "GSN_ClimTier1"] <- gsn[, , "GSN_ClimTier1"] * landMismatch[, , "GSN_ClimTier1"]
    gsn[, , "GSN_ClimTier2"] <- gsn[, , "GSN_ClimTier2"] * landMismatch[, , "GSN_ClimTier2"]

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation and ESA CCI:
    gsn <- toolCorrectOpenEcosystemMismatch(gsn, landUse2020)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LanduseInitialisation data
      totalForest <- dimSums(landUse9[, , c("primforest", "secdforest")], dim = 3) # nolint
      primforestShr <- landUse9[, , "primforest"] / setNames(totalForest + 1e-10, NULL)
      secdforestShr <- landUse9[, , "secdforest"] / setNames(totalForest + 1e-10, NULL)
      # where LanduseInitialisation does not report forest, but we find forest land in
      # GSN data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in GSN data set
      primforest <- setNames(primforestShr, NULL) * gsn[, , paste(getItems(gsn, dim = 3.1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * gsn[, , paste(getItems(gsn, dim = 3.1), "forest", sep = ".")]

      out <- mbind(
        gsn[, , c("crop", "past")],
        new.magpie(getCells(gsn), getYears(gsn), paste(getItems(gsn, dim = 3.1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getItems(gsn, dim = 3.1), "primforest", sep = ".")),
        setNames(secdforest, paste(getItems(gsn, dim = 3.1), "secdforest", sep = ".")),
        new.magpie(getCells(gsn), getYears(gsn), paste(getItems(gsn, dim = 3.1), "urban", sep = "."), fill = 0),
        gsn[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      past <- new.magpie(
        cells_and_regions = getCells(gsn),
        years = getYears(gsn),
        names = paste(getItems(gsn, dim = 3.1), "past", sep = "."),
        fill = 0
      )
      range <- gsn[, , paste(getItems(gsn, dim = 3.1), "past", sep = ".")]

      # separate other land into primary and secondary
      totalOther <- dimSums(landUse9[, , c("primother", "secdother")], dim = 3) # nolint
      primotherShr <- landUse9[, , "primother"] / setNames(totalOther + 1e-10, NULL)
      secdotherShr <- landUse9[, , "secdother"] / setNames(totalOther + 1e-10, NULL)
      # where LanduseInitialisation does not report other land, but we find other land in
      # GSN data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * gsn[, , paste(getItems(gsn, dim = 3.1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * gsn[, , paste(getItems(gsn, dim = 3.1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        past,
        setNames(range, paste(getItems(gsn, dim = 3.1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getItems(gsn, dim = 3.1), "primother", sep = ".")),
        setNames(secdother, paste(getItems(gsn, dim = 3.1), "secdother", sep = "."))
      )
    }
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
