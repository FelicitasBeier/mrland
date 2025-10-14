#' @title calc2ndBioDem
#' @description calculates 2nd generation bioenergy demand
#' @return magpie object with results on country level, weight on country level, unit and description.
#'
#' @param datasource source to be used
#' @param rev data revision the output will be produced for (numeric_version)
#'
#' @examples
#' \dontrun{
#' calcOutput("2ndBioDem")
#' }
#' @import magclass
#' @importFrom madrat readSource calcOutput
#' @importFrom magclass collapseNames time_interpolate mbind lowpass

calc2ndBioDem <- function(datasource, rev = numeric_version("0.1")) {

  if (datasource == "REMIND") {
    x <- readSource("REMIND",
                    subtype = paste0("extensive_",
                                     rev,
                                     "_",
                                     "Primary Energy Production|Biomass|Energy Crops (EJ/yr)"))
    x <- x * 10^3
    x <- collapseNames(x)
    firstRemindYear <- sort(getYears(x))[1]
    x <- time_interpolate(x, seq(1995, 2150, 5), extrapolation_type = "constant")

    # set values in initial years that are not existing in REMIND data to zero
    x[, getYears(x) < firstRemindYear, ] <- 0

    description <- "2nd generation bioenergy demand for different scenarios taken from R2M41 coupled runs"

  } else if (datasource == "Strefler2021") {
    x <- readSource("Strefler2021", subtype = paste0("extensive_", rev))
    x <- x[, , "Primary Energy Production|Biomass|Energy Crops (EJ/yr)"] * 10^3
    x <- collapseNames(x)
    firstRemindYear <- sort(getYears(x))[1]
    x <- time_interpolate(x, seq(1995, 2150, 5), extrapolation_type = "constant")

    # set values in initial years that are not existing in REMIND data to zero
    x[, getYears(x) < firstRemindYear, ] <- 0

    description <- paste("2nd generation bioenergy demand for different scenarios taken",
                         "from Strefler et al 2021 (DOI 10.1038/s41467-021-22211-2)")

  } else if (datasource == "REMMAG") {
    x <- readSource("REMMAG", "biodem")
    # harmonize historic period
    x[, c(1995, 2000, 2005, 2010), ] <- collapseNames(x[, c(1995, 2000, 2005, 2010), "SSP2-Ref-SPA0"])
    description <- "2nd generation bioenergy demand for different scenarios taken from R17M3 coupled runs"

  } else if (datasource == "SSPResults") {
    x <- readSource("SSPResults")
    x <- collapseNames(x[, , "Primary Energy|Biomass|Energy Crops (EJ/yr)"]) * 10^3
    description <- "2nd generation bioenergy demand for different scenarios taken from IIASA SSP database"

  } else if (datasource == "S4N_project") {
    # Total bioenergy demand (including 1st genation, 2nd generation and residues)
    # at country level from IMAGE for 2 different SSP2 scenarios starting in 2005 (in EJ per year)
    imageBe <- readSource("S4Nproject_input", subtype = "bioenergy", convert = "onlycorrect")

    # Transform units: from EJ to PJ
    imageBe <- imageBe * 1e3

    # 1st gen BE demand in MAgPIE in scenario selected for Sim4Nexus (in PJ/yr)
    be1st <- calcOutput("1stBioDem", years = seq(2005, 2100, by = 5), aggregate = FALSE)
    be1st <- collapseNames(be1st[, , "const2030"])
    be1st <- dimSums(be1st, dim = 3)
    # 2nd gen residues in MAgPIE in scenario selected for Sim4Nexus (in PJ/yr)
    res    <- calcOutput("ResFor2ndBioengery", products = "kres", product_aggr = TRUE,
                         add_off = TRUE, years = seq(2005, 2100, by = 5), aggregate = FALSE)
    res    <- collapseNames(res[, , "ssp2"])

    # 2nd generation bioenergy demand: Total BE (IMAGE) - 1st BE (MAgPIE) - residues (MAgPIE)
    imageBe <- imageBe - be1st - res

    # Correct negative values
    imageBe[imageBe < 0] <- 0

    # Fill missing years (1995, 2000) with 2nd generation bioenergy demand from REMMAG data
    remmagBE  <- readSource("REMMAG", subtype = "biodem")
    remmagBE  <- collapseNames(remmagBE[, c("y1995", "y2000"), "SSP2-26-SPA2"])
    x <- new.magpie(getCells(imageBe), paste0("y", seq(1995, 2100, by = 5)), getNames(imageBe))

    x[, getYears(remmagBE), "SSP2"]            <- remmagBE
    x[, getYears(remmagBE), "SSP2_SPA2_26I_D"] <- remmagBE
    x[, getYears(imageBe), ]                   <- imageBe

    # fill missing years in the future
    x <- time_interpolate(x, seq(1995, 2150, 5), extrapolation_type = "constant")

    description <- "2nd generation bioenergy demand for different scenarios provided by IMAGE"

  } else if (datasource == "SSP_and_REM") {
    ssp <- calcOutput("2ndBioDem", datasource = "SSPResults", aggregate = FALSE, rev = rev)
    rem <- calcOutput("2ndBioDem", datasource = "REMIND", aggregate = FALSE, rev = rev)
    if (rev > numeric_version("4.58")) {
      strefler <- calcOutput("2ndBioDem", datasource = "Strefler2021", aggregate = FALSE, rev = rev)
    } else {
      strefler <- NULL
    }

    ssp <- time_interpolate(ssp, getYears(rem), extrapolation_type = "constant")
    x <- mbind(ssp, rem, strefler)

    # years
    years <- getYears(x, as.integer = TRUE)
    yrHist <- years[years > 1995 & years <= 2020]
    yrFut <- years[years >= 2020]

    # apply lowpass filter (not applied on 1st time step, applied separately on historic and future period)
    # only on sceanrios that were existing before revision 4.188
    iter <- 3
    modelsBeforeRev4188 <- "SSPDB|R2M41|R21M42|R32M46|PIK_NPI|PIK_HOS|PIK_HBL|PIK_OPT|PIK_H2C|PIK_GDP|PIK_LIN"
    scenariosBeforeRev4188 <- grep(modelsBeforeRev4188, getItems(x, dim = 3), value = TRUE)

    x[, , scenariosBeforeRev4188] <- mbind(x[, 1995, scenariosBeforeRev4188],
                                           lowpass(x[, yrHist, scenariosBeforeRev4188], i = iter),
                                           lowpass(x[, yrFut, scenariosBeforeRev4188],
                                                   i = iter)[, -1, scenariosBeforeRev4188])

    # sort scenarios alphabetically
    x <- x[, , sort(getNames(x))]

    description <- paste("2nd generation bioenergy demand for different scenarios",
                         "taken from R2M41 coupled runs and from IIASA SSP database")

  } else {
    stop("Unknown datasource", datasource)
  }

  return(list(x = x, weight = NULL,
              description = description,
              unit = "PJ per year",
              note = paste("bioenergy is demanded in the country which is expected",
                           "to produce the bioenergy (demand after trade)")))
}
