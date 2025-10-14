#' @title calcYieldsCalibrated
#' @description This functions calibrates extracted yields from LPJmL to
#'              FAO country level yields
#'
#' @param source        Defines LPJmL version for main crop inputs and isimip replacement.
#'                      For isimip choose crop model/gcm/rcp/co2 combination formatted like this:
#'                      "yields:EPIC-IIASA:ukesm1-0-ll:ssp585:default:3b"
#' @param climatetype   switch between different climate scenarios
#' @param selectyears   Years to be returned (for memory reasons)
#' @param refYear       reference year for calibration
#' @param refYields     assumption for baseline yields with respect to multiple cropping (e.g.,
#'                      FALSE: single-cropped LPJmL yields used as baseline to calculate country-level yields,
#'                      "TRUE:actual:irrig_crop": multicropped yields
#'                                                where LandInG reports current multiple cropping
#'                                                (irrigation- and crop-specific))
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (mask can be:
#'                      "none": no mask applied (only for development purposes)
#'                      "actual:total": currently multicropped areas calculated from total harvested areas
#'                                      and total physical areas per cell from readLanduseLandInG
#'                      "actual:crop" (crop-specific), "actual:irrigation" (irrigation-specific),
#'                      "actual:irrig_crop" (crop- and irrigation-specific),
#'                      "potential:endogenous": potentially multicropped areas given
#'                                              temperature and productivity limits
#'                      "potential:exogenous": potentially multicropped areas given
#'                                             GAEZ suitability classification)
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param areaSource    data source for croparea used in calculation: FAO or LandInG
#'                      Note: when calibrating multicropped yields, LandInG croparea should be used
#' @param average       Averaging period in years used in calcFAOYield (if NULL no averaging is used)
#' @param aggregation   Aggregation at which calibration is performed (options: "country", "continents", "GLO")
#' @param marginal_land Defines which share of marginal land should be included (see options below) and
#'                      whether suitable land under irrigated conditions ("irrigated"),
#'                      under rainfed conditions ("rainfed")
#'                      or suitability under rainfed conditions
#'                      including currently irrigated land (rainfed_and_irrigated)
#'                      should be used. Options combined via ":"
#'                      The different marginal land options are:
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
#' @return magpie object in cellular resolution from reference year onwards
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsCalibrated", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass getYears getNames dimSums mbind
#' @importFrom madrat calcOutput toolConditionalReplace
#' @importFrom mstools toolCoord2Isocell
#' @importFrom withr local_options

calcYieldsCalibrated <- function(datasource = c(lpjml = "ggcmi_phase3_nchecks_9ca735cb", isimip = NULL),
                                 climatetype = "GSWP3-W5E5:historical",
                                 refYear = "y1995", selectyears = seq(1965, 2100, by = 5),
                                 cells = "lpjcell",
                                 multicropping = FALSE, refYields = FALSE,
                                 areaSource = "FAO", marginal_land = "magpie", # nolint
                                 average = 5, aggregation = "country") {

  local_options(magclass_sizeLimit = 1e+12)

  # correct ref year format
  if (!grepl("y", refYear)) {
    refYear <- paste0("y", refYear)
  }

  # extract crop list
  crops          <- setdiff(findset("kcr"), c("betr", "begr"))

  # read FAO and LPJmL yields
  yieldFAOiso    <- calcOutput("FAOYield", cut = 0.98, areaSource = areaSource,
                               average = average,
                               aggregate = FALSE)[, refYear, crops]
  yieldLPJmLgrid <- calcOutput("Yields", datasource = datasource, climatetype = climatetype,
                               selectyears = selectyears,
                               multicropping = multicropping, marginal_land = marginal_land,
                               aggregate = FALSE, supplementary = TRUE)
  yieldLPJmLbase <- calcOutput("Yields", datasource = datasource, climatetype = climatetype,
                               selectyears = selectyears,
                               multicropping = refYields, marginal_land = marginal_land,
                               aggregate = FALSE, supplementary = FALSE)

  years          <- getYears(yieldLPJmLgrid$x, as.integer = TRUE)
  years          <- years[years >= as.integer(gsub("y", "", refYear))]
  weight         <- yieldLPJmLgrid$weight

  otherYields    <- yieldLPJmLgrid$x[, years, setdiff(getItems(yieldLPJmLgrid$x, dim = "crop"), crops)]
  yieldLPJmLgrid <- yieldLPJmLgrid$x[, years, crops]
  yieldLPJmLbase <- yieldLPJmLbase[, refYear, crops]

  # crop-specific cropland area split by irrigation and rainfed
  if (areaSource == "FAO") {

    cropareaMAGgrid <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                                  cellular = TRUE,  cells = "lpjcell",
                                  irrigation = TRUE, aggregate = FALSE)[, refYear, crops]


    # total irrigated & rainfed cropland (for correction of 0 cropland areas)
    proxyMAGgrid    <- dimSums(cropareaMAGgrid, dim = "MAG")

  } else if (areaSource == "LandInG") {

    cropareaMAGgrid <- calcOutput("CropareaLandInG", sectoral = "kcr", physical = TRUE,
                                  irrigation = TRUE, selectyears = refYear,
                                  cellular = TRUE, aggregate = FALSE)[, , crops]
    cropareaMAGgrid <- dimOrder(cropareaMAGgrid, perm = c(2, 1), dim = 3)
    # total irrigated & rainfed cropland (for correction of 0 cropland areas)
    proxyMAGgrid    <- dimSums(cropareaMAGgrid, dim = "crop")

  }

  # Aggregate to country values
  # Crop-specific total cropland area per country
  cropareaMAGiso <- dimSums(cropareaMAGgrid, dim = c("x", "y", "irrigation"))

  # Averaged LPJmL yield per country (LPJmL production / area)
  yieldLPJmLiso  <- dimSums(dimSums(yieldLPJmLbase * cropareaMAGgrid,
                                    dim = 3.2),
                            dim = c("x", "y")) / cropareaMAGiso

  # Correction where no historical crop-specific areas given
  yieldLPJmLiso[cropareaMAGiso == 0] <- (dimSums(dimSums(yieldLPJmLbase * proxyMAGgrid,
                                                         dim = 3.2),
                                                 dim = c("x", "y")) / dimSums(cropareaMAGiso,
                                                                              dim = 3))[cropareaMAGiso == 0]

  # Correction NAs
  yieldLPJmLiso <- toolConditionalReplace(yieldLPJmLiso, "is.na()", 0)

  # Harmonize countries
  yieldLPJmLiso <- yieldLPJmLiso[intersect(getCells(yieldLPJmLiso), getCells(yieldFAOiso)), , ]
  yieldFAOiso   <- yieldFAOiso[intersect(getCells(yieldLPJmLiso), getCells(yieldFAOiso)), , ]

  # Optional: Aggregate to continental or global value
  if (aggregation == "GLO") {
    rel <- data.frame(iso = getItems(yieldFAOiso, dim = "iso"),
                      glo = rep("GLO", length(getItems(yieldFAOiso, dim = "iso"))))
    yieldLPJmLiso <- toolAggregate(yieldLPJmLiso, rel = rel, weight = cropareaMAGiso + 1e-10, from = "iso", to = "glo")
    yieldFAOiso   <- toolAggregate(yieldFAOiso, rel = rel, weight = cropareaMAGiso + 1e-10, from = "iso", to = "glo")
  } else if (grepl("continent", aggregation)) {
    rel <- toolGetMapping("country2continent.csv", where = "mrland")
    rel <- rel[rel$iso %in% getItems(yieldFAOiso, dim = "iso"), ]
    tmp <- toolGetMappingCoord2Country()
    rel <- merge(tmp, rel, by = "iso")
    rel <- rel[order(match(rel$coords, tmp$coords)), ]

    yieldLPJmLiso <- toolAggregate(yieldLPJmLiso, rel = rel, weight = cropareaMAGiso + 1e-10, from = "iso", to = "continent")
    yieldFAOiso   <- toolAggregate(yieldFAOiso, rel = rel, weight = cropareaMAGiso + 1e-10, from = "iso", to = "continent")

    getItems(yieldLPJmLgrid, dim = 1, raw = TRUE) <- paste(rel$coords, rel$iso, rel$continent, sep = ".")
    names(dimnames(yieldLPJmLgrid))[1] <- "x.y.iso.iso1"
  }

  # Yield calibration of LPJmL yields to FAO country yield levels
  out <- toolPatternScaling(yieldLPJmLgrid, yieldLPJmLiso, yieldFAOiso, refYear = refYear)
  # correct dimensions
  if (length(strsplit(names(dimnames(out))[1], split = "\\.")[[1]]) > 3) {
    out <- collapseDim(out, dim = 1.4)
  }

  # Combine with pasture, betr, begr yields that were not calibrated
  getCells(out) <- getCells(otherYields)
  out           <- mbind(out, otherYields)

  return(list(x            = out,
              weight       = weight,
              unit         = "t DM per ha physical area",
              description  = "Calibrated crop yields by plant type and irrigation",
              isocountries = FALSE))
}
