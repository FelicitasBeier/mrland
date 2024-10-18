#' @title calcYields
#'
#' @description This function extracts yields from LPJmL
#'              and transforms them to MAgPIE crops calibrating proxy crops
#'              to FAO yields. Optionally, ISIMIP yields can be returned.
#'
#' @param source        Defines LPJmL version for main crop inputs and isimip replacement.
#'                      For isimip choose crop model/gcm/rcp/co2 combination formatted like this:
#'                      "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b"
#' @param climatetype   Switch between different climate scenarios
#' @param cells         if cellular is TRUE: "magpiecell" for 59199 cells or "lpjcell" for 67420 cells
#' @param selectyears   Years to be returned
#' @param weighting     use of different weights (totalCrop (default),
#'                      totalLUspecific, cropSpecific, crop+irrigSpecific,
#'                      avlCropland, avlCropland+avlPasture)
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (mask can be:
#'                      "none": no mask applied (only for development purposes)
#'                      "actual:total": currently multicropped areas calculated from total harvested areas
#'                                      and total physical areas per cell from readLandInG
#'                      "actual:crop" (crop-specific), "actual:irrigation" (irrigation-specific),
#'                      "actual:irrig_crop" (crop- and irrigation-specific),
#'                      "potential:endogenous": potentially multicropped areas given
#'                                              temperature and productivity limits
#'                      "potential:exogenous": potentially multicropped areas given
#'                                             GAEZ suitability classification)
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param indiaYields   if TRUE returns scaled yields for rainfed crops in India
#' @param scaleFactor   integer value by which indiaYields will be scaled
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
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("Yields", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass getYears add_columns dimSums time_interpolate
#' @importFrom madrat toolFillYears toolGetMapping toolTimeAverage
#' @importFrom mrlandcore toolLPJmLVersion
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom stringr str_split
#' @importFrom withr local_options

calcYields <- function(source = c(lpjml = "ggcmi_phase3_nchecks_9ca735cb", isimip = NULL), # nolint
                       climatetype = "GSWP3-W5E5:historical", cells = "lpjcell",
                       selectyears = seq(1965, 2100, by = 5),
                       weighting = "totalCrop", multicropping = FALSE,
                       indiaYields = FALSE, scaleFactor = 0.3,
                       marginal_land = "magpie") { # nolint

  # Extract argument information
  areaMask      <- paste(str_split(multicropping, ":")[[1]][2],
                         str_split(multicropping, ":")[[1]][3], sep = ":")
  multicropping <- as.logical(str_split(multicropping, ":")[[1]][1])

  # Set up size limit
  local_options(magclass_sizeLimit = 1e+12)

  # LPJmL yields
  yields  <- setYears(calcOutput("YieldsLPJmL", lpjml = source[["lpjml"]], # nolint: undesirable_function_linter.
                                 climatetype = climatetype,
                                 years = selectyears,
                                 cells = cells, aggregate = FALSE),
                      selectyears)
  getSets(yields)["d3.1"] <- "crop"
  getSets(yields)["d3.2"] <- "irrigation"
  # Mapping LPJmL to MAgPIE crops
  lpj2mag   <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")

  if (multicropping) {

    increaseFactor <- calcOutput("MulticroppingYieldIncrease",
                                 lpjml = source[["lpjml"]], # nolint: undesirable_function_linter.
                                 climatetype = climatetype,
                                 selectyears = selectyears,
                                 aggregate = FALSE)

    if (cells == "magpiecell") {
      increaseFactor <- toolCoord2Isocell(increaseFactor)
    }

    # Main-season yield
    mainYield <- yields
    # Off-season yield
    offYield  <- yields * increaseFactor

    # LPJmL to MAgPIE crops
    mainYield <- toolAggregate(mainYield, lpj2mag, dim = 3.1, partrel = TRUE,
                               from = "LPJmL5", to = "MAgPIE")
    offYield  <- toolAggregate(offYield, lpj2mag, dim = 3.1, partrel = TRUE,
                               from = "LPJmL5", to = "MAgPIE")

    # Multiple cropping suitability
    if (areaMask == "none") {
      suitMC <- calcOutput("MulticroppingCells",
                           sectoral = "kcr",
                           scenario = "potential:exogenous",
                           lpjml = c(crop = source[["lpjml"]]), # nolint: undesirable_function_linter.
                           climatetype = climatetype,
                           selectyears = selectyears,
                           aggregate = FALSE)
      # Add pasture to suitMC object with suitability set to 0
      suitMC <- add_columns(suitMC, dim = 3.1, addnm = "pasture", fill = 0)
      # multiple cropping is allowed everywhere
      suitMC[, , ] <- 1
    } else {
      suitMC <- calcOutput("MulticroppingCells", scenario = areaMask,
                           sectoral = "kcr",
                           lpjml =  c(crop = source[["lpjml"]]), # nolint: undesirable_function_linter.
                           climatetype = climatetype,
                           selectyears = selectyears,
                           aggregate = FALSE)
      # Add pasture to suitMC object with suitability set to 0
      suitMC <- add_columns(suitMC, dim = 3.1, addnm = "pasture", fill = 0)
    }

    # Whole year yields under multicropping (main-season yield + off-season yield)
    yields <- mainYield + offYield * suitMC

  } else {

    # LPJmL to MAgPIE crops
    yields  <- toolAggregate(yields, lpj2mag, dim = 3.1, partrel = TRUE,
                             from = "LPJmL5", to = "MAgPIE")
  }

  # Check for NAs
  if (any(is.na(yields))) {
    stop("produced NA yields")
  }

  # Use FAO data to scale proxy crops to reasonable levels (global, static factor)
  prodFAO  <- collapseNames(calcOutput("FAOmassbalance_pre", aggregate = FALSE)[, , "production"][, , "dm"])
  areaMAg  <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE, aggregate = FALSE)
  faoyears <- intersect(getYears(prodFAO), getYears(areaMAg))

  cropsMAg <- findset("kcr")
  missing  <- c("betr", "begr")
  cropsMAg <- setdiff(cropsMAg, missing)
  prodFAO  <- add_columns(prodFAO[, , cropsMAg], addnm = missing, dim = 3.1)
  prodFAO[, , missing] <- 0

  yieldsFAO <- dimSums(prodFAO[, faoyears, ], dim = 1) / dimSums(areaMAg[, faoyears, ], dim = 1)
  yieldsFAO <- setYears(toolTimeAverage(yieldsFAO[, 1993:1997, ], 5))
  calib     <- new.magpie(cells_and_regions = "GLO",
                          years = NULL,
                          names = c(getNames(yieldsFAO), "pasture"),
                          fill = 1,
                          sets = c("iso", "year", "data"))
  calib[, , "oilpalm"]   <- yieldsFAO[, , "oilpalm"] / yieldsFAO[, , "groundnut"]   # LPJmL proxy is groundnut
  calib[, , "cottn_pro"] <- yieldsFAO[, , "cottn_pro"] / yieldsFAO[, , "groundnut"] # LPJmL proxy is groundnut
  calib[, , "foddr"]     <- yieldsFAO[, , "foddr"] / yieldsFAO[, , "maiz"]          # LPJmL proxy is maize
  calib[, , "others"]    <- yieldsFAO[, , "others"] / yieldsFAO[, , "maiz"]         # LPJmL proxy is maize
  calib[, , "potato"]    <- yieldsFAO[, , "potato"] / yieldsFAO[, , "sugr_beet"]    # LPJmL proxy is sugarbeet

  # Recalibrate yields for proxies
  yields <- yields * calib[, , getNames(yields, dim = 1)]

  if (!is.na(source["isimip"])) { # nolint: undesirable_function_linter.
    isimipYields <- calcOutput("ISIMIP3bYields", subtype = source[["isimip"]], # nolint: undesirable_function_linter.
                               cells = cells, aggregate = FALSE)
    commonVars  <- intersect(getNames(yields), getNames(isimipYields))
    commonYears <- intersect(getYears(yields), getYears(isimipYields))

    #  harmonize to LPJml
    cfg       <- toolLPJmLVersion(version = source["lpjml"], # nolint: undesirable_function_linter.
                                  climatetype = climatetype)
    repHarmon <- toolHarmonize2Baseline(x = isimipYields[, commonYears, commonVars],
                                        base = yields[, commonYears, commonVars],
                                        ref_year = cfg$ref_year_gcm)
    gc()
    # convert to array for memory
    yields    <- as.array(yields)
    repHarmon <- as.array(repHarmon)
    yields[, commonYears, commonVars] <- repHarmon[, commonYears, commonVars]
    yields    <- as.magpie(yields, spatial = 1)
    repHarmon <- as.magpie(repHarmon, spatial = 1)

  }

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

      cropAreaWeight <- new.magpie(cells_and_regions = getCells(yields),
                                   years = NULL,
                                   names = getNames(yields),
                                   fill = NA)
      cropAreaWeight[, , findset("kcr")] <- crop + 10e-10
      cropAreaWeight[, , "pasture"]      <- mbind(setNames(past + 10e-10, "irrigated"),
                                                  setNames(past + 10e-10, "rainfed"))

    } else if (weighting == "cropSpecific") {

      cropAreaWeight <- new.magpie(cells_and_regions = getCells(yields),
                                   years = NULL,
                                   names = getNames(yields, dim = 1),
                                   fill = NA)

      cropAreaWeight[, , findset("kcr")] <- dimSums(crop, dim = 3.1) + 10e-10
      cropAreaWeight[, , "pasture"]      <- past + 10e-10

    } else {

      cropAreaWeight <- new.magpie(cells_and_regions = getCells(yields),
                                   years = NULL,
                                   names = getNames(yields, dim = 1),
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

    cropAreaWeight <- new.magpie(cells_and_regions = getCells(yields),
                                 years = NULL,
                                 names = getNames(yields, dim = 1),
                                 fill = avlCrop)

    cropAreaWeight[, , "pasture"] <- pmax(avlCrop,
                                          dimSums(lu1995[, , c("primforest", "secdforest", "forestry", "past")],
                                                  dim = 3)) + 10e-10

  } else {

    stop("Weighting setting is not available.")
  }

  if (any(is.na(cropAreaWeight))) stop("NAs in weights.")

  # Special case for India case study
  if (indiaYields) {
    yields["IND", , "rainfed"] <- yields["IND", , "rainfed"] * scaleFactor
  }

  return(list(x            = yields,
              weight       = cropAreaWeight,
              unit         = "t per ha",
              description  = "Yields in tons per hectar for different crop types.",
              isocountries = FALSE))
}
