#' @title calcAvlCropland
#'
#' @description Calculates the total available cropland per grid cell,
#'              based on physical cropland suitability data or other criteria,
#'              such as constraints on cropland expansion
#'
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
#'                        "no_marginal:rainfed_and_irrigated"
#'                        in a magclass object to be used as magpie input.
#' }
#' @param cell_upper_bound Upper bound for cropland at the grid cell level.
#'                         Even if, for instance, the total available cropland area equals the land area in a grid cell,
#'                         cropland cannot be expanded above this value.
#' @param cells            (deprecated) only lpjcell (67420 cells)
#' @param country_level    Whether output shall be at country level.
#'                         Requires aggregate=FALSE in calcOutput.
#' @param luhBaseYear      Base year of LUH land area
#'
#' @return magpie object in cellular resolution
#' @author Patrick v. Jeetze, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("AvlCropland", aggregate = FALSE)
#' }
#'
#' @importFrom madrat readSource calcOutput
#' @importFrom magclass dimSums getCells getYears getNames mbind collapseDim as.magpie
#' @importFrom magpiesets addLocation
#' @importFrom mstools toolCoord2Isocell toolAggregateCell2Country
#'

calcAvlCropland <- function(marginal_land = "magpie", cell_upper_bound = 0.9, country_level = FALSE, # nolint
                            cells = "lpjcell", luhBaseYear = "y1995") {

  # extract function arguments
  marginalLand <- marginal_land # nolint
  cellUpperBound <- cell_upper_bound # nolint
  countryLevel <- country_level # nolint

  if (is.numeric(luhBaseYear)) {
    luhBaseYear <- paste0("y", luhBaseYear)
  }

  # read luh data in chosen base year
  luh <- calcOutput("LUH3", landuseTypes = "magpie", aggregate = FALSE,
                    cellular = TRUE, irrigation = FALSE, years = luhBaseYear)
  # sum land area per grid cell
  landarea <- dimSums(luh, dim = 3)

  x <- as.magpie(NULL)

  if (any(grepl("all_marginal", marginalLand)) || marginalLand == "magpie") {

    if (marginalLand == "magpie") {
      cropsuit <- readSource("Zabel2014", subtype = paste("all_marginal", "rainfed_and_irrigated",
                                                          sep = ":"),
                             convert = "onlycorrect")
    } else {
      cropsuit <- readSource("Zabel2014", subtype = paste("all_marginal",
                                                          unlist(strsplit(marginalLand, split = ":"))[2],
                                                          sep = ":"),
                             convert = "onlycorrect")
    }

    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "all_marginal"
    x <- mbind(x, tmp)

  }

  if (any(grepl("q33_marginal", marginalLand)) || marginalLand == "magpie") {


    if (marginalLand == "magpie") {
      cropsuit <- readSource("Zabel2014", subtype = paste("q33_marginal", "rainfed_and_irrigated",
                                                          sep = ":"),
                             convert = "onlycorrect")
    } else {
      cropsuit <- readSource("Zabel2014", subtype = paste("q33_marginal",
                                                          unlist(strsplit(marginalLand, split = ":"))[2],
                                                          sep = ":"),
                             convert = "onlycorrect")
    }
    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "q33_marginal"
    x <- mbind(x, tmp)

  }

  if (any(grepl("q50_marginal", marginalLand))) {

    cropsuit <- readSource("Zabel2014", subtype = paste("q50_marginal",
                                                        unlist(strsplit(marginalLand, split = ":"))[2],
                                                        sep = ":"),
                           convert = "onlycorrect")
    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "q50_marginal"
    x <- mbind(x, tmp)

  }


  if (any(grepl("q66_marginal", marginalLand))) {

    cropsuit <- readSource("Zabel2014", subtype = paste("q66_marginal",
                                                        unlist(strsplit(marginalLand, split = ":"))[2],
                                                        sep = ":"),
                           convert = "onlycorrect")

    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "q66_marginal"
    x <- mbind(x, tmp)

  }

  if (any(grepl("q75_marginal", marginalLand))) {

    cropsuit <- readSource("Zabel2014", subtype = paste("q75_marginal",
                                                        unlist(strsplit(marginalLand, split = ":"))[2],
                                                        sep = ":"),
                           convert = "onlycorrect")

    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "q75_marginal"
    x <- mbind(x, tmp)

  }

  if (any(grepl("no_marginal", marginalLand)) || marginalLand == "magpie") {

    if (marginalLand == "magpie") {
      cropsuit <- readSource("Zabel2014", subtype = paste("no_marginal", "rainfed_and_irrigated",
                                                          sep = ":"),
                             convert = "onlycorrect")
    } else {
      cropsuit <- readSource("Zabel2014", subtype = paste("no_marginal",
                                                          unlist(strsplit(marginalLand, split = ":"))[2],
                                                          sep = ":"),
                             convert = "onlycorrect")
    }
    # make sure that suitable cropland is not larger than total land area
    cropsuit <- pmin(cropsuit, landarea)
    # set upper bound for cropland at grid cell level
    # in each grid cell cropland cannot be expanded above this threshold
    cropsuit <- cropsuit * cellUpperBound
    # cropland suitability is corrected where LUH reports (more) cropland
    cropsuit <- pmax(cropsuit, luh[, , "crop"])

    tmp <- cropsuit
    getNames(tmp) <- "no_marginal"
    x <- mbind(x, tmp)

  }

  if (countryLevel) {

    out <- toolCountryFill(dimSums(x, dim = c("x", "y")), fill = 0)

  } else {

    if (cells == "lpjcell") {

      out <- x

    } else {
      stop("Please specify cells argument")
    }
  }

  return(list(x            = out,
              weight       = NULL,
              unit         = "Mha",
              description  = "Cropland suitability based on Zabel et al. (2014)
                              with different suitability thresholds
                              ('all_marginal', 'q33_marginal', 'no_marginal').",
              isocountries = FALSE))
}
