#' @title calcOzoneYieldShock
#' @description calculate Ozone yield shocks
#' Data from the EAT-Lancet deepdive on Ozone shock effects on crop yields.
#' @param weighting     use of different weights (totalCrop (default),
#'                      totalLUspecific, cropSpecific, crop+irrigSpecific,
#'                      avlCropland, avlCropland+avlPasture)
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
#' @author Jake Tommey
#' @examples
#' \dontrun{
#' calcOutput("OzoneYieldShock")
#' }
#' @importFrom mstools toolGetMappingCoord2Country

calcOzoneYieldShock <- function(
  weighting = "totalCrop",
  marginal_land = "magpie" #nolint
) {
  shocks <- readSource("OzoneYieldShock", convert = "onlycorrect")
  shocks <- add_columns(shocks, addnm = "y2020", dim = 2, fill = 0)
  # time interpolate the shock
  shocks <- magclass::time_interpolate(
    shocks,
    interpolated_year = seq(1965, 2150, 5),
    integrate_interpolated_years = TRUE,
    extrapolation_type = "constant"
  )

  mapping <- toolGetMappingCoord2Country()
  mapping$coordiso <- paste(mapping$coords, mapping$iso, sep = ".")
  shocks <- toolAggregate(shocks, mapping, from = "iso", to = "coordiso")
  getSets(shocks) <- c("x.y.iso", "year", "rcp.crop")

  # create a MAgPIE object with the correct shape.
  cropNames <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")$MAgPIE
  yieldShock <- new.magpie(
    cells_and_regions = getCells(shocks),
    years = getYears(shocks),
    names = c(paste0("rcp2p6.", cropNames), paste0("rcp7p0.", cropNames))
  )

  yieldShock[, , "tece"] <- shocks[, , "Wheat"]
  yieldShock[, , "soybean"] <- shocks[, , "soybean"]
  yieldShock[, , "rice_pro"] <- shocks[, , "Rice"]
  yieldShock[, , "sugr_beet"] <- shocks[, , "sugarbeet"]
  yieldShock[, , "maiz"] <- shocks[, , "Maize"]
  otherCrops <- yieldShock[, , c("tece", "soybean", "rice_pro", "sugr_beet", "maiz"), invert = TRUE]
  yieldShock[, , getItems(otherCrops, dim = 3)] <- (dimSums(shocks, dim = 3) / length(getItems(shocks, dim = 3)))
  yieldShock[, , "pasture"] <- 0

  getSets(yieldShock) <- c("x", "y", "iso", "year", "rcp", "crop")

  cropAreaWeights <- calcOutput(
    "YieldsWeight",
    weighting = weighting,
    marginal_land = marginal_land,
    aggregate = FALSE
  )

  return(
    list(
      x = yieldShock,
      weight = cropAreaWeights,
      unit = "%",
      description = "percentage yield shock due to ozone",
      isocountries = FALSE
    )
  )
}
