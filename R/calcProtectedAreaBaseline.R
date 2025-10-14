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
#' \code{\link{calcProtectArea}}
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
    luh3 <- calcOutput("LUH3",
      landuseTypes = "LUH3", aggregate = FALSE,
      cellular = TRUE,
      yrs = c(1995, 2000, 2005, 2010, 2015, 2020)
    )

    # extend the data set to all time steps provided in the protected area data
    # i.e. use the data from the year 2015 for the year 2020.
    luh3 <- setCells(
      luh3,
      getCells(PABaseline)
    )

    # calculate total land area
    landArea <- dimSums(luh3[, "y1995", ], dim = 3)

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

    if (nclasses %in% c("seven", "nine")) {
      # differentiate primary and secondary forest based on luh3 data
      totforestluh <- dimSums(luh3[, , c("primf", "secdf")], dim = 3)
      primforestShr <- luh3[, , "primf"] / setNames(totforestluh + 1e-10, NULL)
      secdforestShr <- luh3[, , "secdf"] / setNames(totforestluh + 1e-10, NULL)
      # where luh2 does not report forest, but we find forest land in
      # protected area data, set share of secondary forest land to 1
      secdforestShr[secdforestShr == 0 & primforestShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with
      # land pools in protected area data set
      primforest <- setNames(primforestShr, NULL) * PABaseline[, , "forest"]
      secdforest <- setNames(secdforestShr, NULL) * PABaseline[, , "forest"]
      primforest <- setNames(primforest, sub("forest", "primforest", getNames(primforest)))
      secdforest <- setNames(secdforest, sub("forest", "secdforest", getNames(secdforest)))

      out <- mbind(
        PABaseline[, , c("crop", "past")],
        new.magpie(
          cells_and_regions = getCells(PABaseline),
          years = getYears(PABaseline),
          names = paste0(c("WDPA", "WDPA_I-II-III", "WDPA_IV-V-VI"), ".forestry"),
          fill = 0
        ),
        primforest,
        secdforest,
        new.magpie(
          cells_and_regions = getCells(PABaseline),
          years = getYears(PABaseline),
          names = paste0(c("WDPA", "WDPA_I-II-III", "WDPA_IV-V-VI"), ".urban"),
          fill = 0
        ),
        PABaseline[, , "other"]
      )
    } else {
      stop("Option specified for argument 'nclasses' does not exist.")
    }

    if (nclasses == "nine") {
      # separate pasture into pasture and rangeland
      totgrassluh <- dimSums(luh3[, , c("pastr", "range")], dim = 3)
      pastShr <- luh3[, , "pastr"] / setNames(totgrassluh + 1e-10, NULL)
      rangeShr <- luh3[, , "range"] / setNames(totgrassluh + 1e-10, NULL)
      # where luh2 does not report grassland, but we find grassland in
      # protected area data, set share of rangeland to 1
      rangeShr[pastShr == 0 & rangeShr == 0] <- 1
      # multiply shares of pasture and rangeland with pasture in protected area data
      past <- setNames(pastShr, NULL) * PABaseline[, , "past"]
      range <- setNames(rangeShr, NULL) * PABaseline[, , "past"]
      range <- setNames(range, sub("past", "range", getNames(range)))

      # separate other land into primary and secondary
      tototherluh <- dimSums(luh3[, , c("primn", "secdn")], dim = 3)
      primotherShr <- luh3[, , "primn"] / setNames(tototherluh + 1e-10, NULL)
      secdotherShr <- luh3[, , "secdn"] / setNames(tototherluh + 1e-10, NULL)
      # where luh2 does not report other land, but we find other land in
      # protected area data, set share of secondary other land to 1
      secdotherShr[secdotherShr == 0 & primotherShr == 0] <- 1
      # multiply shares of primary and secondary non-forest veg with other land
      primother <- setNames(primotherShr, NULL) * PABaseline[, , "other"]
      secdother <- setNames(secdotherShr, NULL) * PABaseline[, , "other"]
      primother <- setNames(primother, sub("other", "primother", getNames(primother)))
      secdother <- setNames(secdother, sub("other", "secdother", getNames(secdother)))

      out <- mbind(
        out[, , "crop"],
        past,
        range,
        out[, , c("forestry", "primforest", "secdforest", "urban")],
        primother,
        secdother
      )
    }
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
