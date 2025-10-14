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
    luh3 <- calcOutput("LUH3",
      landuseTypes = "LUH3", aggregate = FALSE,
      cellular = TRUE,
      yrs = 2015
    )
    getYears(luh3) <- NULL
    getCells(luh3) <- getCells(ic)

    # calculate total land area
    landArea <- dimSums(luh3, dim = 3)

    # urban land
    urbanLand <- calcOutput("UrbanLandFuture",
      subtype = "LUH3", aggregate = FALSE,
      timestep = "5year", cells = "lpjcell"
    )

    # make sure that irrecoverable carbon land is not greater than total land area minus urban area
    landNoUrban <- setYears(landArea, "y2020") - setCells(urbanLand[, "y2020", "SSP2"], getCells(landArea))
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

    if (nclasses %in% c("seven", "nine")) {

      # differentiate primary and secondary forest based on luh3 data
      totForestLUH <- dimSums(luh3[, , c("primf", "secdf")], dim = 3) # nolint
      primforestShr <- luh3[, , "primf"] / setNames(totForestLUH + 1e-10, NULL)
      secdforestShr <- luh3[, , "secdf"] / setNames(totForestLUH + 1e-10, NULL)
      # where luh2 does not report forest, but we find forest land in
      # irrecoverable carbon data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in irrecoverable carbon data set
      primforest <- setNames(primforestShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "forest", sep = ".")]

      out <- mbind(
        ic[, , c("crop", "past")],
        new.magpie(getCells(ic), getYears(ic), paste(getNames(ic, dim = 1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getNames(ic, dim = 1), "primforest", sep = ".")),
        setNames(secdforest, paste(getNames(ic, dim = 1), "secdforest", sep = ".")),
        new.magpie(getCells(ic), getYears(ic), paste(getNames(ic, dim = 1), "urban", sep = "."), fill = 0),
        ic[, , "other"]
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
      # irrecoverable carbon data, set share of rangeland to 1
      rangeShr[pastShr == 0 & rangeShr == 0] <- 1
      # multiply shares of pasture and rangeland with pasture in irrecoverable carbon data
      past <- setNames(pastShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "past", sep = ".")]
      range <- setNames(rangeShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "past", sep = ".")]

      # separate other land into primary and secondary
      totOtherLUH <- dimSums(luh3[, , c("primn", "secdn")], dim = 3) # nolint
      primotherShr <- luh3[, , "primn"] / setNames(totOtherLUH + 1e-10, NULL)
      secdotherShr <- luh3[, , "secdn"] / setNames(totOtherLUH + 1e-10, NULL)
      # where luh2 does not report other land, but we find other land in
      # irrecoverable carbon data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * ic[, , paste(getNames(ic, dim = 1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        setNames(past, paste(getNames(ic, dim = 1), "past", sep = ".")),
        setNames(range, paste(getNames(ic, dim = 1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getNames(ic, dim = 1), "primother", sep = ".")),
        setNames(secdother, paste(getNames(ic, dim = 1), "secdother", sep = "."))
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
