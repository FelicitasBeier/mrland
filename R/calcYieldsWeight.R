#' @title calcYieldsWeight
#'
#' @description This function calculates the crop area weightings to use for yields.
#' @param weighting     use of different weights (totalCrop (default),
#'                      totalLUspecific, cropSpecific, crop+irrigSpecific,
#'                      avlCropland, avlCropland+avlPasture)
#' @param cells         if cellular is TRUE: "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
#' @param marginal_land  Defines which share of marginal land should be included (see options below) and
#'                whether suitable land under irrigated conditions ("irrigated"), under rainfed conditions ("rainfed")
#'                or suitability under rainfed conditions including currently irrigated land (rainfed_and_irrigated)
#'                should be used. Options combined via ":"
#'                The different marginal land options are:
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
#' @importFrom mrlandcore toolLPJmLVersion
#' @importFrom stringr str_split
#' @importFrom withr local_options

calcYieldsWeight <- function(cells = "lpjcell", weighting = "totalCrop", marginal_land = "magpie") { # nolint

  yieldNames <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")$MAgPIE
  isos <- toolGetMappingCoord2Country()
  yieldCells <- paste(isos$coords, isos$iso, sep = ".")


  # Weight for spatial aggregation
  if (weighting == "totalCrop") {

    cropAreaWeight <- dimSums(calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = FALSE,
                                         cellular = TRUE, cells = cells, aggregate = FALSE,
                                         years = "y1995", round = 6),
                              dim = 3) + 10e-10

  } else if (weighting %in% c("totalLUspecific", "cropSpecific", "crop+irrigSpecific")) {

    crop <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE, irrigation = TRUE,
                       cellular = TRUE, cells = cells, aggregate = FALSE, years = "y1995", round = 6)

    past <- calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                       input_magpie = TRUE, cells = cells, years = "y1995", round = 6)[, , "past"]

    if (weighting == "crop+irrigSpecific") {

      cropAreaWeight <- new.magpie(cells_and_regions = yieldCells,
                                   years = NULL,
                                   names = c(paste(yieldNames, "irrigated", sep = "."),
                                             paste(yieldNames, "rainfed", sep = ".")),
                                   fill = NA)
      cropAreaWeight[, , findset("kcr")] <- crop + 10e-10
      cropAreaWeight[, , "pasture"]      <- mbind(setNames(past + 10e-10, "irrigated"),
                                                  setNames(past + 10e-10, "rainfed"))

    } else if (weighting == "cropSpecific") {

      cropAreaWeight <- new.magpie(cells_and_regions = yieldCells,
                                   years = NULL,
                                   names = yieldNames,
                                   fill = NA)

      cropAreaWeight[, , findset("kcr")] <- dimSums(crop, dim = 3.1) + 10e-10
      cropAreaWeight[, , "pasture"]      <- past + 10e-10

    } else {

      cropAreaWeight <- new.magpie(cells_and_regions = yieldCells,
                                   years = NULL,
                                   names = yieldNames,
                                   fill = (dimSums(crop, dim = 3) + 10e-10))

      cropAreaWeight[, , "pasture"] <- past + 10e-10

    }

  } else if (weighting == "avlCropland") {

    cropAreaWeight <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land, cells = cells,
                                          country_level = FALSE, aggregate = FALSE),
                               NULL) + 10e-10

  } else if (weighting == "avlCropland+avlPasture") {

    avlCrop <- setNames(calcOutput("AvlCropland", marginal_land = marginal_land, cells = cells,
                                   country_level = FALSE, aggregate = FALSE), "avlCrop")

    lu1995  <- setYears(calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, nclasses = "seven",
                                   input_magpie = TRUE, cells = cells, years = "y1995", round = 6), NULL)

    cropAreaWeight <- new.magpie(cells_and_regions = yieldCells,
                                 years = NULL,
                                 names = yieldNames,
                                 fill = avlCrop)

    cropAreaWeight[, , "pasture"] <- pmax(avlCrop,
                                          dimSums(lu1995[, , c("primforest", "secdforest", "forestry", "past")],
                                                  dim = 3)) + 10e-10

  } else {

    stop("Weighting setting is not available.")
  }

  if (any(is.na(cropAreaWeight))) stop("NAs in weights.")

  return(list(x            = cropAreaWeight,
              weight       = NULL,
              unit         = "Mha",
              description  = "Yields in tons per hectar for different crop types.",
              isocountries = FALSE))
}
