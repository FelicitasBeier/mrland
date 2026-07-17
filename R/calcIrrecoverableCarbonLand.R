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
#' \code{\link{readNoon2022}}
#'
#' @examples
#' \dontrun{
#' calcOutput("calcIrrecoverableCarbonLand", aggregate = FALSE)
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
    .alignLuDims <- function(lu, aData) {
      getYears(lu) <- NULL
      getCells(lu) <- getCells(aData)
      return(lu)
    }

    baseLandUse <- calcOutput("LanduseInitialisation", nclasses = "five", cellular = TRUE,
                              input_magpie = TRUE, aggregate = FALSE, years = 2020)
    landUse2020 <- .alignLuDims(baseLandUse[, 2020, ], ic)

    # calculate total land area
    landArea <- dimSums(landUse2020, dim = 3)

    # urban land from LanduseInitialisation
    landUse9 <- calcOutput("LanduseInitialisation",
      nclasses = "nine", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE,
      years = 2020
    )[, "y2020", ]
    landUse9 <- .alignLuDims(landUse9, ic)
    urbanLand2020 <- landUse9[, , "urban"]

    # make sure that irrecoverable carbon land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand2020
    getYears(landNoUrban) <- getYears(ic)
    # compute mismatch factor
    icTotalLand <- mbind(
      dimSums(ic[, , "IrrC_50pc"], dim = 3.2),
      dimSums(ic[, , "IrrC_75pc"], dim = 3.2),
      dimSums(ic[, , "IrrC_95pc"], dim = 3.2),
      dimSums(ic[, , "IrrC_99pc"], dim = 3.2),
      dimSums(ic[, , "IrrC_75pc_30by30"], dim = 3.2),
      dimSums(ic[, , "IrrC_95pc_30by30"], dim = 3.2),
      dimSums(ic[, , "IrrC_99pc_30by30"], dim = 3.2)
    )
    landMismatch <- setNames(landNoUrban, NULL) / icTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
    # correct irrecoverable carbon data
    ic[, , "IrrC_50pc"] <- ic[, , "IrrC_50pc"] * landMismatch[, , "IrrC_50pc"]
    ic[, , "IrrC_75pc"] <- ic[, , "IrrC_75pc"] * landMismatch[, , "IrrC_75pc"]
    ic[, , "IrrC_95pc"] <- ic[, , "IrrC_95pc"] * landMismatch[, , "IrrC_95pc"]
    ic[, , "IrrC_99pc"] <- ic[, , "IrrC_99pc"] * landMismatch[, , "IrrC_99pc"]
    ic[, , "IrrC_75pc_30by30"] <- ic[, , "IrrC_75pc_30by30"] * landMismatch[, , "IrrC_75pc_30by30"]
    ic[, , "IrrC_95pc_30by30"] <- ic[, , "IrrC_95pc_30by30"] * landMismatch[, , "IrrC_95pc_30by30"]
    ic[, , "IrrC_99pc_30by30"] <- ic[, , "IrrC_99pc_30by30"] * landMismatch[, , "IrrC_99pc_30by30"]

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation and ESA CCI:
    ic <- toolCorrectOpenEcosystemMismatch(ic, landUse2020)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LanduseInitialisation data
      totalForest <- dimSums(landUse9[, , c("primforest", "secdforest")], dim = 3) # nolint
      primforestShr <- landUse9[, , "primforest"] / setNames(totalForest + 1e-10, NULL)
      secdforestShr <- landUse9[, , "secdforest"] / setNames(totalForest + 1e-10, NULL)
      # where LanduseInitialisation does not report forest, but we find forest land in
      # irrecoverable carbon data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in irrecoverable carbon data set
      primforest <- setNames(primforestShr, NULL) * ic[, , paste(getItems(ic, dim = 3.1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * ic[, , paste(getItems(ic, dim = 3.1), "forest", sep = ".")]

      out <- mbind(
        ic[, , c("crop", "past")],
        new.magpie(getCells(ic), getYears(ic), paste(getItems(ic, dim = 3.1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getItems(ic, dim = 3.1), "primforest", sep = ".")),
        setNames(secdforest, paste(getItems(ic, dim = 3.1), "secdforest", sep = ".")),
        new.magpie(getCells(ic), getYears(ic), paste(getItems(ic, dim = 3.1), "urban", sep = "."), fill = 0),
        ic[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      past <- new.magpie(
        cells_and_regions = getCells(ic),
        years = getYears(ic),
        names = paste(getItems(ic, dim = 3.1), "past", sep = "."),
        fill = 0
      )
      range <- ic[, , paste(getItems(ic, dim = 3.1), "past", sep = ".")]

      # separate other land into primary and secondary
      totalOther <- dimSums(landUse9[, , c("primother", "secdother")], dim = 3) # nolint
      primotherShr <- landUse9[, , "primother"] / setNames(totalOther + 1e-10, NULL)
      secdotherShr <- landUse9[, , "secdother"] / setNames(totalOther + 1e-10, NULL)
      # where LanduseInitialisation does not report other land, but we find other land in
      # irrecoverable carbon data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * ic[, , paste(getItems(ic, dim = 3.1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * ic[, , paste(getItems(ic, dim = 3.1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        past,
        setNames(range, paste(getItems(ic, dim = 3.1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getItems(ic, dim = 3.1), "primother", sep = ".")),
        setNames(secdother, paste(getItems(ic, dim = 3.1), "secdother", sep = "."))
      )
    }
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
