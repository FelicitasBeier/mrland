#' @title calcYieldsMAgPIEcrops
#'
#' @description
#' Transforms LPJmL yields to MAgPIE crop types (kcr) and applies proxy crop
#' calibration. Optionally applies FAO yield calibration when a calibration
#' list is provided.
#'
#' @param lpjml         LPJmL version for main crop inputs
#' @param climatetype   Climate scenario or historical baseline
#' @param selectyears   Years to be returned
#' @param multicropping Multicropping activated (TRUE) or not (FALSE) and
#'                      Multiple Cropping Suitability mask selected
#'                      (e.g. TRUE:actual:total; TRUE:none; FALSE)
#' @param calibration   NULL for uncalibrated yields (default)
#'                      or a list with FAO calibration settings.
#'                      Expected keys:
#'                      refYear, refYields, areaSource, average, aggregation.
#'                      Example:
#'                      list(refYear = "y1995", refYields = FALSE,
#'                           areaSource = "FAO", average = 5,
#'                           aggregation = "country")
#'
#' @return magpie object in cellular resolution
#'
#' @author Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("YieldsMAgPIEcrops", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass getYears getCells getItems add_columns dimSums
#'   new.magpie collapseDim dimOrder mbind
#' @importFrom madrat toolGetMapping toolTimeAverage toolAggregate
#'   toolConditionalReplace toolPatternScaling
#' @importFrom mstools toolGetMappingCoord2Country
#' @importFrom stringr str_split
#' @importFrom withr local_options

calcYieldsMAgPIEcrops <- function(lpjml = "ggcmi_phase3_nchecks_9ca735cb",
                                  climatetype = "GSWP3-W5E5:historical",
                                  selectyears = seq(1965, 2100, by = 5),
                                  multicropping = FALSE,
                                  calibration = NULL) {

  # Increase object size limit
  local_options(magclass_sizeLimit = 1e+12)

  # LPJmL yields
  yieldLPJmLgrid <- calcOutput("YieldsLPJmL", lpjml = lpjml,
                               climatetype = climatetype,
                               selectyears = selectyears,
                               multicropping = multicropping,
                               supplementary = TRUE,
                               aggregate = FALSE)
  yields  <- setYears(yieldLPJmLgrid$x, selectyears)

  # LPJmL to MAgPIE crops
  lpj2mag <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")
  yields  <- toolAggregate(yields, lpj2mag, dim = "crop", partrel = TRUE,
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
    lpj2mag <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mrlandcore")
    if (lpj2mag$LPJmL5[lpj2mag$MAgPIE == "oilpalm"] == "oil crops groundnut") {
      # To Do: yields[, , "oilpalm"] <- yields1st[, , "oilpalm"] + yields2nd[, , "oilpalm"]
    } else {
      warning("The LPJmL-to-MAgPIE crop mapping has changed and the special treatment
              of oilpalm yields is no longer necessary.
              Please remove this chunk of code from calcYieldsMAgPIEcrops.")
    }
  }

  # Check for NAs
  if (any(is.na(yields))) {
    stop("produced NA yields")
  }

  # Use FAO data to scale proxy crops to reasonable levels (global, static factor)
  #### To Do (discuss with Kristine, Mike, Feli): replace with calcFAOYield,
  #### but needs adjustments in calcFAOYield to allow to return global average FAO yield
  #### rather than regional yield
  #### This is possible already, right? Just use calcOutput(FAOYield, aggregate = "GLO")
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

  description <- "Yields for different MAgPIE crop types"

  ################################
  ### Optional FAO calibration ###
  ################################
  # This section mimicks the yield calibrtion to FAO levels that takes place
  # in MAgPIE's preloop. Different settings can be chosen via argument list.
  if (!is.null(calibration)) {

    # Extract calibration arguments
    refYear     <- calibration$refYear
    refYields   <- calibration$refYields
    areaSource  <- calibration$areaSource
    average     <- calibration$average
    aggregation <- calibration$aggregation

    if (is.null(refYear) | is.null(refYields) | is.null(areaSource)) {
      stop("Please provide a list with: refYear, refYields, areaSource, (optional: average), (optional: aggregation) for calibration of yields.")
    }

    if (!grepl("y", refYear)) {
      refYear <- paste0("y", refYear)
    }

    # extract crop list
    crops <- setdiff(findset("kcr"), c("betr", "begr"))

    # read FAO and LPJmL yields
    yieldFAOiso    <- calcOutput("FAOYield", cut = 0.98, areaSource = areaSource,
                                 average = average,
                                 aggregate = FALSE)[, refYear, crops]

    yieldLPJmLgrid <- yields
    # LPJmL yields for calibration baseline
    yieldLPJmLbase <- setYears(calcOutput("YieldsLPJmL", lpjml = lpjml,
                                          climatetype = climatetype,
                                          selectyears = selectyears,
                                          multicropping = refYields,
                                          supplementary = FALSE,
                                          aggregate = FALSE),
                               selectyears)
    yieldLPJmLbase <- toolAggregate(yieldLPJmLbase, lpj2mag, dim = "crop", partrel = TRUE,
                                    from = "LPJmL5", to = "MAgPIE")

    years          <- getYears(yieldLPJmLgrid, as.integer = TRUE)
    years          <- years[years >= as.integer(gsub("y", "", refYear))]

    otherYields    <- yieldLPJmLgrid[, years, setdiff(getItems(yieldLPJmLgrid, dim = "crop"), crops)]
    yieldLPJmLgrid <- yieldLPJmLgrid[, years, crops]
    yieldLPJmLbase <- yieldLPJmLbase[, refYear, crops]

    # crop-specific cropland area split by irrigation and rainfed
    if (areaSource == "FAO") {
      cropareaMAGgrid <- calcOutput("Croparea", sectoral = "kcr", physical = TRUE,
                                    cellular = TRUE, cells = "lpjcell",
                                    irrigation = TRUE, aggregate = FALSE)[, refYear, crops]
      proxyMAGgrid    <- dimSums(cropareaMAGgrid, dim = "MAG")
    } else if (areaSource == "LandInG") {
      cropareaMAGgrid <- calcOutput("CropareaLandInG", sectoral = "kcr", physical = TRUE,
                                    irrigation = TRUE, selectyears = refYear,
                                    cellular = TRUE, aggregate = FALSE)[, , crops]
      cropareaMAGgrid <- dimOrder(cropareaMAGgrid, perm = c(2, 1), dim = 3)
      proxyMAGgrid    <- dimSums(cropareaMAGgrid, dim = "crop")
    } else {
      stop("Please specify which area should be used for calculation.\n         Note: LandInG should be FAO-consistent.")
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
    yields        <- mbind(out, otherYields)
    description   <- "Yields (calibrated to FAO country level values) for different MAgPIE crop types"
  }

  return(list(x            = yields,
              weight       = NULL,  ### Kristine, Jan: Should we provide a weight here?
              unit         = "tDM per ha",
              description  = description,
              isocountries = FALSE))
}
