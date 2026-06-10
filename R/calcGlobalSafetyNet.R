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
    luh3 <- calcOutput("LUH3",
      landuseTypes = "LUH3", aggregate = FALSE,
      cellular = TRUE, irrigation = FALSE,
      yrs = 2015
    )
    getYears(luh3) <- NULL
    getCells(luh3) <- getCells(gsn)

    # calculate total land area
    landArea <- dimSums(luh3, dim = 3)

    # urban land
    urbanLand <- calcOutput("UrbanLandFuture",
      subtype = "LUH3", aggregate = FALSE,
      timestep = "5year", cells = "lpjcell"
    )
    urbanLand <- setCells(urbanLand, getCells(landArea))

    # make sure that GSN land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand[, "y2015", "SSP2"]
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
    # initialisation (LUH) and ESA CCI:
    luIni <- calcOutput("LanduseInitialisation",
      nclasses = "seven", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE
    )[, "y2020", ]
    getYears(luIni) <- NULL
    getCells(luIni) <- getCells(gsn)
    gsn <- toolCorrectOpenEcosystemMismatch(gsn, luIni)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LUH3 data
      totForestLUH <- dimSums(luh3[, , c("primf", "secdf")], dim = 3) # nolint
      primforestShr <- luh3[, , "primf"] / setNames(totForestLUH + 1e-10, NULL)
      secdforestShr <- luh3[, , "secdf"] / setNames(totForestLUH + 1e-10, NULL)
      # where luh3 does not report forest, but we find forest land in
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
      totOtherLUH <- dimSums(luh3[, , c("primn", "secdn")], dim = 3) # nolint
      primotherShr <- luh3[, , "primn"] / setNames(totOtherLUH + 1e-10, NULL)
      secdotherShr <- luh3[, , "secdn"] / setNames(totOtherLUH + 1e-10, NULL)
      # where luh3 does not report other land, but we find other land in
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
