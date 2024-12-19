#' @title calcYields
#'
#' @description This function extracts yields from calcYieldsLPJmL
#'              and transforms them to MAgPIE crops calibrating proxy crops
#'              to FAO yields. Optionally, ISIMIP yields can be returned.
#'
#' @param datasource    Choose LPJmL version for main crop inputs and optionally ISIMIP version
#'                      For ISIMIP choose crop model/gcm/rcp/co2 combination formatted like this:
#'                      "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b"
#' @param climatetype   Switch between different climate scenarios
#' @param selectyears   Years to be returned
#' @param weighting     Use of different weights (totalCrop (default),
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
#'                                              growing conditions from LPJmL
#'                      "potential:exogenous": potentially multicropped areas given
#'                                             GAEZ suitability classification)
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param indiaYields    If TRUE returns scaled yields for rainfed crops in India
#' @param scaleFactor    Integer value by which indiaYields will be scaled
#' @param marginal_land  Defines which share of marginal land should be included (see options below) and
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
#' @importFrom mrlandcore toolLPJmLHarmonize
#' @importFrom mstools toolHarmonize2Baseline
#' @importFrom stringr str_split
#' @importFrom withr local_options

calcYields <- function(datasource = c(lpjml = "ggcmi_phase3_nchecks_9ca735cb", isimip = NULL),
                       climatetype = "GSWP3-W5E5:historical",
                       selectyears = seq(1965, 2100, by = 5),
                       weighting = "totalCrop", multicropping = FALSE,
                       indiaYields = FALSE, scaleFactor = 0.3,
                       marginal_land = "magpie") { # nolint

  # Set up size limit
  local_options(magclass_sizeLimit = 1e+12)

  # LPJmL yields
  yields  <- setYears(calcOutput("YieldsLPJmL", lpjml = datasource[["lpjml"]],
                                 climatetype = climatetype,
                                 selectyears = selectyears,
                                 multicropping = multicropping,
                                 aggregate = FALSE),
                      selectyears)

  # LPJmL to MAgPIE crops
  lpj2mag   <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")
  yields  <- toolAggregate(yields, lpj2mag, dim = 3.2, partrel = TRUE,
                           from = "LPJmL5", to = "MAgPIE")

  # Perennial crops in MAgPIE
  if (as.logical(str_split(multicropping, ":")[[1]][1])) {
    ### To Do (discuss with Kristine and Jens): Do we want a special treatment for oilpalm.
    ### If so (To Do: Feli): adjust calcYieldsLPJmL such that groundnut first and second
    ### season yield for groundnut can be returned and used here.
    ### Also (To Do): This chunk of code is actually valid
    ### not only for the multiple cropping case, but generally.
    ### Once the multiple cropping functionality is fully tested, the if-condition
    ### should be removed.

    # The MAgPIE perennial crop "oilpalm" is grown throughout the whole year
    # but proxied with an LPJmL crop with seasonality ("groundnut").
    # Therefore, the main season yield should be adjusted to the whole year yield
    if (lpj2mag$LPJmL5[lpj2mag$MAgPIE == "oilpalm"] == "oil crops groundnut") {
      # To Do: yields[, , "oilpalm"] <- yields1st[, , "oilpalm"] + yields2nd[, , "oilpalm"]
    } else {
      warning("The LPJmL-to-MAgPIE crop mapping has changed and the special treatment
              of oilpalm yields is no longer necessary.
              Please remove this chunk of code from calcYields.")
    }
  }

  # Check for NAs
  if (any(is.na(yields))) {
    stop("produced NA yields")
  }

  # Use FAO data to scale proxy crops to reasonable levels (global, static factor)
  #### To Do (Kristine, Mike, Feli): replace with calcFAOYield,
  #### but needs adjustments in calcFAOYield to allow to return global average FAO yield
  #### rather than regional yield
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
  ### Question (Kristine): Do we need special treatment for oilpalm here as well?
  ### Due to fact that it is a perennial in MAgPIE proxied by groundnut from LPJmL
  calib[, , "cottn_pro"] <- yieldsFAO[, , "cottn_pro"] / yieldsFAO[, , "groundnut"] # LPJmL proxy is groundnut
  calib[, , "foddr"]     <- yieldsFAO[, , "foddr"] / yieldsFAO[, , "maiz"]          # LPJmL proxy is maize
  calib[, , "others"]    <- yieldsFAO[, , "others"] / yieldsFAO[, , "maiz"]         # LPJmL proxy is maize
  calib[, , "potato"]    <- yieldsFAO[, , "potato"] / yieldsFAO[, , "sugr_beet"]    # LPJmL proxy is sugarbeet

  # Recalibrate yields for proxies
  yields <- yields * calib[, , getItems(yields, dim = "crop")]

  if (!is.na(datasource["isimip"])) {
    isimipYields <- calcOutput("ISIMIP3bYields", subtype = datasource[["isimip"]],
                               cells = "lpjcell", aggregate = FALSE)
    commonVars  <- intersect(getNames(yields), getNames(isimipYields))
    commonYears <- intersect(getYears(yields), getYears(isimipYields))

    #  harmonize to LPJml
    cfg       <- toolLPJmLHarmonize(lpjmlversion = datasource["lpjml"],
                                    climatetype  = climatetype)
    repHarmon <- toolHarmonize2Baseline(x = isimipYields[, commonYears, commonVars],
                                        base = yields[, commonYears, commonVars],
                                        ref_year = cfg$refYearGcm)
    gc()
    # convert to array for memory
    yields    <- as.array(yields)
    repHarmon <- as.array(repHarmon)
    yields[, commonYears, commonVars] <- repHarmon[, commonYears, commonVars]
    yields    <- as.magpie(yields, spatial = 1)
    repHarmon <- as.magpie(repHarmon, spatial = 1)

    #### Question (Jens, Kristine, Jan): Would we want to be able to apply
    #### the multiple cropping logic on ISIMIP yields as well?
    #### If so: I would suggest to use the same yield increase factor
    #### (as calculated from LPJmL grass GPP)
    #### and move the multicropping yield increase calculation from calcYieldsLPJmL
    #### to a separate tool-/calc-Function

  }

  # select weight for crop and pasture yields
  cropAreaWeight <- calcOutput("YieldsWeight",
                               weighting = weighting,
                               marginal_land = marginal_land,
                               aggregate = FALSE)

  # Special case for India case study
  if (indiaYields) {
    yields["IND", , "rainfed"] <- yields["IND", , "rainfed"] * scaleFactor
    #### Question (Kristine, Misko, Edna): Should this scale factor be applied to both
    #### LPJmL and ISIMIP yields? I guess Vartika based this on LPJmL yields only?
    #### Also: should we check whether this is needed with the new LPJmL data still?
    #### And finally: I guess it won't be necessary with the multiple cropping update anymore.
    ####              Should we deprecate this argument once multiple cropping is merged?
  }

  return(list(x            = yields,
              weight       = cropAreaWeight,
              unit         = "tDM per ha",
              description  = "Yields for different MAgPIE crop types",
              isocountries = FALSE))
}
