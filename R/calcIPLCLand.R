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
#' \code{\link{readLandMark}}
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
    .alignLuDims <- function(lu, aData) {
      getYears(lu) <- NULL
      getCells(lu) <- getCells(aData)
      return(lu)
    }

    baseLandUse <- calcOutput("LanduseInitialisation", nclasses = "five", cellular = TRUE,
                              input_magpie = TRUE, aggregate = FALSE, years = 2020)
    landUse2020 <- setYears(.alignLuDims(baseLandUse, iplc), "y2020")

    # calculate total land area
    landArea <- dimSums(landUse2020, dim = 3)

    # urban land from LanduseInitialisation
    landUse9 <- setYears(calcOutput("LanduseInitialisation",
      nclasses = "nine", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE,
      years = 2020
    ), "y2020")
    landUse9 <- .alignLuDims(landUse9, iplc)
    urbanLand2020 <- landUse9[, , "urban"]

    # make sure that IPLC land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand2020
    getYears(landNoUrban) <- getYears(iplc)

    # compute mismatch factor
    iplcTotalLand <- dimSums(iplc, dim = 3.2)
    landMismatch <- setNames(landNoUrban, NULL) / iplcTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)

    # correct IPLC data
    iplc <- iplc * landMismatch

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation and ESA CCI:
    iplc <- toolCorrectOpenEcosystemMismatch(iplc, landUse2020)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LanduseInitialisation data
      totalForest <- dimSums(landUse9[, , c("primforest", "secdforest")], dim = 3) # nolint
      primforestShr <- landUse9[, , "primforest"] / setNames(totalForest + 1e-10, NULL)
      secdforestShr <- landUse9[, , "secdforest"] / setNames(totalForest + 1e-10, NULL)
      # where LanduseInitialisation does not report forest, but we find forest land in
      # IPLC data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary forest with
      # land pools in IPLC data set
      primforest <- setNames(primforestShr, NULL) * iplc[, , paste(getItems(iplc, dim = 3.1), "forest", sep = ".")]
      secdforest <- setNames(secdforestShr, NULL) * iplc[, , paste(getItems(iplc, dim = 3.1), "forest", sep = ".")]

      out <- mbind(
        iplc[, , c("crop", "past")],
        new.magpie(getCells(iplc), getYears(iplc), paste(getItems(iplc, dim = 3.1), "forestry", sep = "."), fill = 0),
        setNames(primforest, paste(getItems(iplc, dim = 3.1), "primforest", sep = ".")),
        setNames(secdforest, paste(getItems(iplc, dim = 3.1), "secdforest", sep = ".")),
        new.magpie(getCells(iplc), getYears(iplc), paste(getItems(iplc, dim = 3.1), "urban", sep = "."), fill = 0),
        iplc[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      past <- new.magpie(
        cells_and_regions = getCells(iplc),
        years = getYears(iplc),
        names = paste(getItems(iplc, dim = 3.1), "past", sep = "."),
        fill = 0
      )
      range <- iplc[, , paste(getItems(iplc, dim = 3.1), "past", sep = ".")]

      # separate other land into primary and secondary
      totalOther <- dimSums(landUse9[, , c("primother", "secdother")], dim = 3) # nolint
      primotherShr <- landUse9[, , "primother"] / setNames(totalOther + 1e-10, NULL)
      secdotherShr <- landUse9[, , "secdother"] / setNames(totalOther + 1e-10, NULL)
      # where LanduseInitialisation does not report other land, but we find other land in
      # IPLC data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * iplc[, , paste(getItems(iplc, dim = 3.1), "other", sep = ".")]
      secdother <- setNames(secdotherShr, NULL) * iplc[, , paste(getItems(iplc, dim = 3.1), "other", sep = ".")]

      out <- mbind(
        out[, , "crop"],
        past,
        setNames(range, paste(getItems(iplc, dim = 3.1), "range", sep = ".")),
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        setNames(primother, paste(getItems(iplc, dim = 3.1), "primother", sep = ".")),
        setNames(secdother, paste(getItems(iplc, dim = 3.1), "secdother", sep = "."))
      )
    }
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
