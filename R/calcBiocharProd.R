#' @title calcBiocharProd
#' @description calculates biochar production
#'
#' @param datasource source to be used (either "REMIND" or "Dorndorf2026")
#' @param rev data revision the output will be produced for (numeric_version).
#' @return magpie object with results on country level, weight on country level, unit and description.
#'
#' @author Isabelle Weindl
#' @seealso
#' \code{\link{calc2ndBioDem}}, \code{\link{calcGHGPrices}}
#' @examples
#' \dontrun{
#' calcOutput("BiocharProd")
#' }
#' @importFrom madrat readSource
#' @importFrom magclass collapseNames time_interpolate mbind complete_magpie

calcBiocharProd <- function(datasource, rev = numeric_version("0.1")) {

  if (datasource == "REMIND") {
  
    indicator <- "SE|Biochar"
    unitInfo  <- "PJ per year"
    unit      <- "EJ/yr"

    # Mapping between MAgPIE-internal names for pyrolysis plant set-ups and
    # REMIND technology keys based on the history of REMIND reporting names,
    # where the full REMIND variable name is: indicator|<key> (unit).
    techMaps <- c(biopyrchp  = "w/ heat and power",
                  biopyrhe   = "w/ heat",
                  biopyrliq  = "w/ liquids",
                  biopyronly = "w/o co-product")

    remindSubtype <- function(key) {
      paste0("extensive_", rev, "_", paste0(indicator, key, " (", unit, ")"))
    }

    x <- do.call(mbind, lapply(names(techMaps), function(tech) {
      tmp <- readSource("REMIND", subtype = remindSubtype(paste0("|", techMap[[tech]])))
      add_dimension(tmp, dim = 3.1, nm = tech)
    }))
    # sort scenarios alphabetically
    x <- x[, , sort(getNames(x, dim = 2))]

    # unit conversion from EJ to PJ
    x <- x * 10^3
    
    # if data is not available from REMIND reports, set biochar production to zero
    x[is.na(x)] <- 0
    
    x <- collapseNames(x)
    firstRemindYear <- sort(getYears(x))[1]
    x <- time_interpolate(x, seq(1995, 2150, 5), extrapolation_type = "constant")

    # set values in initial years that are not existing in REMIND data to zero
    x[, getYears(x) < firstRemindYear, ] <- 0

    x <- complete_magpie(x, fill = 0)

    description <- paste("biochar production for different scenarios",
                         "taken from REMIND standalone and REMIND-MAgPIE coupled runs")

    
  } else if (datasource == "Dorndorf2026") {

    unitInfo <- "PJ per year"

    # biochar production per plant setup from the published REMIND runs of
    # Dorndorf et al 2026, returned by the read function as bc_plant63.scenario
    x <- readSource("Dorndorf2026", subtype = "Baseline")

    # electricity-only co-product plant is added to the CHP plant
    elec <- x[, , "biopyrelec"]
    getItems(elec, dim = 3.1) <- "biopyrchp"
    x <- mbind(x[, , "biopyrchp"] + elec,
               x[, , c("biopyrhe", "biopyrliq", "biopyronly")])

    # unit conversion from EJ to PJ
    x <- x * 10^3

    # if data is not available from the reports, set biochar production to zero
    x[is.na(x)] <- 0

    # note: no collapseNames here, it would drop the single-element scenario
    # level and its header row in the cs3 file
    firstRemindYear <- sort(getYears(x))[1]
    x <- time_interpolate(x, seq(1995, 2150, 5), extrapolation_type = "constant")

    # set values in initial years that are not existing in the data to zero
    x[, getYears(x) < firstRemindYear, ] <- 0

    x <- complete_magpie(x, fill = 0)

    description <- paste("biochar production for the baseline scenario",
                         "taken from the REMIND runs of Dorndorf et al 2026",
                         "(PB650-BCpess-CTS01-noFuel-CHP17-EJ100)")

  } else if (datasource == "REMIND_and_Dorndorf") {

    unitInfo <- "PJ per year"

    rem <- calcOutput("BiocharProd", datasource = "REMIND",
                      aggregate = FALSE, rev = rev)
    dor <- calcOutput("BiocharProd", datasource = "Dorndorf2026",
                      aggregate = FALSE, rev = rev)

    dor <- time_interpolate(dor, getYears(rem), extrapolation_type = "constant")

    x <- mbind(rem, dor)
    # sort scenarios alphabetically
    x <- x[, , sort(getNames(x))]

    description <- paste("biochar production for different scenarios taken from",
                         "REMIND standalone and REMIND-MAgPIE coupled runs and",
                         "from the REMIND runs of Dorndorf et al 2026")

  } else {
    stop("Unknown datasource", datasource)
  }

  return(list(
    x = x,
    weight = NULL,
    description = description,
    unit = unitInfo
  ))
}
