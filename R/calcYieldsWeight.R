#' @title calcYieldsWeight
#'
#' @description This function calculates the crop area weightings to use for yields.
#'
#' @param weighting     area-based weighting options for yield aggregation
#'                      ("totalCrop" (default),
#'                      "totalLUspecific", "cropSpecific", "crop+irrigSpecific",
#'                      "avlCropland", "avlCropland+potentiallyIrrigatedAreas",
#'                      "avlCropland+avlPasture")
#' @param selectyears   Years to be returned
#' @param lpjml         lpjml version, only required if potentially irrigated areas
#'                      are used as aggregation weight
#' @param climatetype   different climate scenarios, only required if potentially
#'                      irrigated areas are used as aggregation weight
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (mask can be:
#'                      "none": no mask applied (only for development purposes)
#'                      "actual:total": currently multicropped areas calculated from total harvested areas
#'                                      and total physical areas per cell from readLandInG
#'                      "actual:crop" (crop-specific), "actual:irrigation" (irrigation-specific),
#'                      "actual:irrig_crop" (crop- and irrigation-specific),
#'                      "potential:endogenous": potentially multicropped areas given
#'                                              growing conditions from LPJmL
#'                      "potential:exogenous": potentially multicropped areas given
#'                                             GAEZ suitability classification)
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param marginal_land  Defines which share of marginal land (see options below)
#'                       should be included if avlCropland is chosen as weighting option and
#'                       whether suitable land under irrigated conditions ("irrigated"),
#'                       under rainfed conditions ("rainfed")
#'                       or suitability under rainfed conditions including
#'                       currently irrigated land ("rainfed_and_irrigated")
#'                       should be used. Options combined via ":"
#'                       The different marginal land options are:
#' \itemize{
#' \item \code{"all_marginal"}: All marginal land (suitability index between 0-0.33) is included as suitable
#' \item \code{"q33_marginal"}: The bottom tertile (suitability index below 0.13) of the
#' marginal land area is excluded.
#' \item \code{"q50_marginal"}: The bottom  half (suitability index below 0.18) of the
#' marginal land area is excluded.
#' \item \code{"q66_marginal"}: The first and second tertile (suitability index below 0.23) of
#' the marginal land area are excluded.
#' \item \code{"q75_marginal"}: The first, second and third quartiles (suitability index below 0.25)
#' of the marginal land are are excluded
#' \item \code{"no_marginal"}: Areas with a suitability index of 0.33 and lower are excluded.
#' \item \code{"magpie"}: Returns "all_marginal:rainfed_and_irrigated",
#'                        "q33_marginal:rainfed_and_irrigated" and
#'                        "no_marginal:rainfed_and_irrigated" in a magclass object to be used as magpie input.
#' }
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsWeight", yields, aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass getYears add_columns dimSums time_interpolate
#' @importFrom madrat toolFillYears toolGetMapping toolTimeAverage
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom stringr str_split
#' @importFrom withr local_options

calcYieldsWeight <- function(weighting = "totalCrop", lpjml = NULL, climatetype = NULL,
                             multicropping = NULL, selectyears = seq(1995, 2100, by = 5),
                             marginal_land = "q33_marginal:rainfed_and_irrigated") { # nolint

  # extract dimension information
  yieldNames <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral",
                               where = "mappingfolder")$MAgPIE
  isos       <- toolGetMappingCoord2Country()
  yieldCells <- paste(isos$coords, isos$iso, sep = ".")

  # magpie crop types
  kcr <- magpiesets::findset("kcr")

  # object with correct dimensions
  cropAreaWeight <- add_dimension(new.magpie(cells_and_regions = yieldCells,
                                             years = NULL,
                                             names = yieldNames,
                                             fill = NA),
                                  dim = 3.2,
                                  add = "irrigation",
                                  nm = c("rainfed", "irrigated"))

  #########################################################
  ############ Weight for spatial aggregation #############
  #########################################################
  #### Current cropland as basis for aggregation ####
  if (weighting == "totalCrop") {

    # croparea in initialization year as weight for all irrigated/rainfed crops and pasture areas
    totalCroparea <- dimSums(calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = FALSE,
                                         cellular = TRUE, aggregate = FALSE,
                                         years = "y1995", round = 6),
                              dim = 3) + 10e-10
    cropAreaWeight[, , ] <- totalCroparea

  } else if (weighting %in% c("totalLUspecific", "cropSpecific", "crop+irrigSpecific")) {
    # crop area in initialization year
    crop <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = TRUE,
                       cellular = TRUE, aggregate = FALSE, years = "y1995", round = 6)
    # pasture area in initialization year
    past <- calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                       input_magpie = TRUE, years = "y1995", round = 6)[, , "past"]

    if (weighting == "crop+irrigSpecific") {

      # every irrigated/rainfed crop is weighted with its specific crop area in the initialization year
      cropAreaWeight[, , kcr] <- crop + 10e-10
      cropAreaWeight[, , "pasture"]      <- mbind(setNames(past + 10e-10, "irrigated"),
                                                  setNames(past + 10e-10, "rainfed"))

    } else if (weighting == "cropSpecific") {

      # every crop is weighted with its specific crop area in the initialization year
      cropAreaWeight[, , kcr] <- dimSums(crop, dim = "irrigation") + 10e-10
      cropAreaWeight[, , "pasture"] <- past + 10e-10

    } else {

      # total croparea as weight for crops
      cropAreaWeight[, , kcr] <- dimSums(crop, dim = 3) + 10e-10
      # pasture area as weight for pasture
      cropAreaWeight[, , "pasture"] <- past + 10e-10

    }

    #### Available cropland as basis for aggregation ####
  } else if (weighting == "avlCropland") {

    # available cropland as weight for all irrigated/rainfed crops and pasture areas
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                    country_level = FALSE, aggregate = FALSE),
                        NULL) + 10e-10

    cropAreaWeight[, , ] <- avlCrop

  } else if (weighting == "avlCropland+potentiallyIrrigatedAreas") {

    # available cropland as weight for all rainfed crops and pasture areas
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                   country_level = FALSE, aggregate = FALSE),
                        NULL) + 10e-10
    # potentially irrigated areas as weight for irrigated crops and pasture areas
    pia <- calcOutput("PotIrrigAreas", cropAggregation = TRUE,
                      lpjml = lpjml, climatetype = climatetype,
                      selectyears = selectyears, iniyear = 1995,
                      multicropping = multicropping,
    # standard options (Question: How to hand them over more elegantly?)
                      efrMethod = "VMF:fair", irrigationsystem = "initialization",
                      accessibilityrule = "CV:2", rankmethod = "USD_m3:GLO:TRUE",
                      gainthreshold = 10, allocationrule = "optimization",
                      yieldcalib = FALSE, comAg = TRUE,
                      fossilGW = TRUE, transDist = 100,
                      landScen = "potCropland:NULL", cropmix = "hist_total",
                      aggregate = FALSE)[, "y1995", "off"][, , "ssp2"]

    cropAreaWeight[, , "rainfed"] <- avlCrop
    cropAreaWeight[, , "irrigated"] <- pia

    # Question: How to handle pasture area? As below? Or treat same as crops?

  } else if (weighting == "avlCropland+avlPasture") {

    # available cropland as weight for all irrigated/rainfed crops
    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land,
                                   country_level = FALSE, aggregate = FALSE),
                        "avlCrop") + 10e-10
    # land use classes in 1995
    lu1995  <- setYears(calcOutput("LanduseInitialisation", aggregate = FALSE,
                                   cellular = TRUE, nclasses = "seven",
                                   input_magpie = TRUE, years = "y1995", round = 6),
                        NULL)

    cropAreaWeight[, , kcr] <- avlCrop

    # weight for pasture is available cropland or other land use classes in 1995 if these exceed available cropland
    cropAreaWeight[, , "pasture"] <- pmax(avlCrop,
                                          dimSums(lu1995[, , c("primforest", "secdforest", "forestry", "past")],
                                                  dim = 3)) + 10e-10

  } else {
    warning("No area-based weight for yield aggregation selected or ",
            "weighting setting is not available. ",
            "Using uniform weight.")
    cropAreaWeight <- add_dimension(new.magpie(cells_and_regions = yieldCells,
                                               years = NULL,
                                               names = yieldNames,
                                               fill = 1),
                                    dim = 3.2,
                                    add = "irrigation",
                                    nm = c("rainfed", "irrigated"))
  }

  if (any(is.na(cropAreaWeight))) stop("NAs in weights.")

  return(list(x            = cropAreaWeight,
              weight       = NULL,
              unit         = "Mha",
              description  = "Area-based weight for yields",
              isocountries = FALSE))
}
