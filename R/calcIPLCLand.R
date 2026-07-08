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
    luh3 <- calcOutput("LUH3",
      landuseTypes = "LUH3", aggregate = FALSE,
      cellular = TRUE, irrigation = FALSE,
      yrs = 2015
    )
    getYears(luh3) <- NULL
    getCells(luh3) <- getCells(iplc)

    # calculate total land area
    landArea <- dimSums(luh3, dim = 3)

    # urban land
    urbanLand <- calcOutput("UrbanLandFuture",
      subtype = "LUH3", aggregate = FALSE,
      timestep = "5year", cells = "lpjcell"
    )
    urbanLand <- setCells(urbanLand, getCells(landArea))

    # make sure that IPLC land is not greater than total land area minus urban area
    landNoUrban <- landArea - urbanLand[, "y2015", "SSP2"]
    getYears(landNoUrban) <- getYears(iplc)

    # compute mismatch factor
    iplcTotalLand <- dimSums(iplc, dim = 3.2)
    landMismatch <- setNames(landNoUrban, NULL) / iplcTotalLand
    landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)

    # correct IPLC data
    iplc <- iplc * landMismatch

    # Consider mismatches in the classification of open
    # ecosystems into pasture and other between land-use
    # initialisation (LUH) and ESA CCI:
    luIni <- calcOutput("LanduseInitialisation",
      nclasses = "seven", aggregate = FALSE, cellular = TRUE, input_magpie = TRUE
    )[, "y2020", ]
    getYears(luIni) <- NULL
    getCells(luIni) <- getCells(iplc)
    iplc <- toolCorrectOpenEcosystemMismatch(iplc, luIni)

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on LUH3 data
      totForestLUH <- dimSums(luh3[, , c("primf", "secdf")], dim = 3) # nolint
      primforestShr <- luh3[, , "primf"] / setNames(totForestLUH + 1e-10, NULL)
      secdforestShr <- luh3[, , "secdf"] / setNames(totForestLUH + 1e-10, NULL)
      # where luh3 does not report forest, but we find forest land in
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
      totOtherLUH <- dimSums(luh3[, , c("primn", "secdn")], dim = 3) # nolint
      primotherShr <- luh3[, , "primn"] / setNames(totOtherLUH + 1e-10, NULL)
      secdotherShr <- luh3[, , "secdn"] / setNames(totOtherLUH + 1e-10, NULL)
      # where luh3 does not report other land, but we find other land in
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
