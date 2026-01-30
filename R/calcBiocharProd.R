#' @title calcBiocharProd
#' @description calculates biochar production
#'
#' @param datasource source to be used
#' @param rev data revision the output will be produced for (numeric_version).
#'
#' @param biocharsystems selects the biochar production systems for which data is returned.
#' Options are disagg which returns data for all available biochar production systems separately,
#' and agg which returns the total of all biochar production systems.
#' @param unit (available: EJ/yr and Mt/yr)
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

calcBiocharProd <- function(datasource, rev = numeric_version("0.1"),
                            biocharsystems = "disagg", unit = "EJ/yr") {

  if (datasource == "REMIND") {

    if (unit == "EJ/yr") {
      indicatorName <- "SE|Biochar"
    } else if (unit == "Mt/yr") {
      indicatorName <- "SE|Biochar Mt"
    }

    if (biocharsystems == "agg") {
      x <- readSource("REMIND",
        subtype = paste0("extensive_",
                         rev,
                         "_",
                         paste0(indicatorName, " (", unit, ")"))
      )

      # sort scenarios alphabetically
      x <- x[, , sort(getNames(x, dim = 1))]

    } else if (biocharsystems == "disagg") {

      # Mapping between MAgPIE-internal names for pyrolysis plant set-ups and
      # REMIND technology keys based on the history of REMIND reporting names,
      # where the full REMIND variable name is: indicatorName|<key> (unit).
      # Note: If future data revisions include REMIND simulations with a different
      # set of biochar technologies, the mapping has to be adapted accordingly.

      oldMap <- c(
        # MAgPIE names =  old REMIND reporting key
        kontiki    = "KonTiki",
        biopyrchp  = "biopyrCHP",
        biopyrelec = "biopyrElec",
        biopyrhe   = "biopyrHeat",
        biopyrliq  = "biopyrFuel",
        biopyronly = "biopyrOnly"
      )

      newMap <- c(
        # MAgPIE names =  new REMIND reporting key
        kontiki    = "kontiki_NA", # not available in REMIND anymore, filled with Zeros
        biopyrchp  = "w/ heat and power",
        biopyrelec = "biopyrelec_NA",  # not available in REMIND anymore, filled with Zeros
        biopyrhe   = "w/ heat",
        biopyrliq  = "w/ liquids",
        biopyronly = "w/o co-product"
      )

      # Version threshold that triggers switch to new REMIND reporting names
      # Current threshold is artificial and needs to be specified once
      # new exogenous scenarios are available.
      revNewStyle <- numeric_version("11.11")

      if (rev >= revNewStyle) {
        specMap <- newMap
      } else {
        specMap <- oldMap
      }

      x <- NULL
      for (magTech in names(specMap)) {
        remTech <- specMap[[magTech]]
        tmp <- readSource("REMIND",
          subtype = paste0("extensive_",
                           rev,
                           "_",
                           paste0(indicatorName, "|", remTech, " (", unit, ")"))
        )
        tmp <- add_dimension(tmp, dim = 3.1, nm = magTech)
        x <- mbind(x, tmp)
      }

      # sort scenarios alphabetically
      x <- x[, , sort(getNames(x, dim = 2))]

    } else {
      stop("Unknown biocharsystems setting", biocharsystems)
    }

    if (unit == "EJ/yr") {
      # unit conversion from EJ to PJ
      x <- x * 10^3
      unitInfo <- "PJ per year"
    } else if (unit == "Mt/yr") {
      unitInfo <- "Mt per year"
    }

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

  } else {
    stop("Unknown datasource", datasource)
  }

  return(list(
    x = x,
    weight = NULL,
    description = description,
    unit = unitInfo,
    note = paste("Biochar production is provided without disaggregation by biomass feedstock type.")
  ))
}
