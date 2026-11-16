#' @title calcBiocharGatePrice
#' @description Calculates time-dependent, conversion-only levelized biochar gate
#' prices for biochar production systems that exclude feedstcok costs. Gate prices
#' are derived from raw REMIND biochar cost parametrization (e.g. CAPEX proxy per
#' annual biomass feedstock throughput, fixed O&M share, variable O&M per tDM input)
#' and converted to USD per GJ biochar using biochar conversion efficiencies and
#' biochar energy content.
#'
#' @param bcSys Selects biochar production systems.
#' Options: "default", or subset of bc_tech63 (e.g. "SP-500", "P-Woolf"). Subsets
#' should only be selected if outputSys is set to "bc_sys".
#' @param variant Cost variants to capture uncertainty range.
#' Options: "low", "central", "high".
#' @param baseYear Base year for learning curves.
#' @param learningRates Named list of learning rates per variant.
#' For Kontiki, no learning is assumed. If set to NULL, defaults are used. If set
#' to FALSE, not learning will be assumed.
#' @param capRecF Capital recovery factor (annualization rate).
#' Default is 0.05 (5%/yr).
#' @param markupHigh Markup (fraction) applied to central gate prices to define high
#' cost variant (if no explicit high data exists). Default is 0.3 (= +30%).
#' @param fallbackMultLow Multiplier applied to inco0_main if inco0_low is missing and
#' variant == "low". Default is 0.7.
#' @param sysCostScaling Multipliers to scale from SP-500 costs to other temperature systems
#' if these are not covered by the raw input data. They should be provided as named
#' vector, e.g. c("SP-400"=0.85,"SP-650"=1.15,"SP-800"=1.35,"FP-500"=2.6,"P-KonTiki"=3.6).
#' If NULL, only systems that can be constructed from input data directly are returned.
#' @param outputSys "bc_sys" (default) returns gate prices aggregated to bc_tech63.
#' "biopyr" returns gate prices for individual REMIND plant set-ups.
#' @return List of magpie objects with results on global level, empty weight,
#' unit and description.
#'
#' @author Isabelle Weindl
#' @seealso [calcBiocharEfficiency()], [calcBiocharAttributes()]
#' @examples
#' \dontrun{
#'   calcOutput("BiocharGatePrice", bcSys = "default")
#' }
calcBiocharGatePrice <- function(
    bcSys = "default",
    variant = c("low", "central", "high"),
    baseYear = "y2025",
    learningRates = NULL,
    capRecF = 0.05,
    markupHigh = 0.30,
    fallbackMultLow = 0.70,
    sysCostScaling =
        c("SP-400" = 0.85, "SP-650" = 1.15, "SP-800" = 1.35,
          "FP-500" = 2.60, "P-KonTiki" = 3.60), outputSys = "bc_sys") {

  ## Defaults for bcSys selection
  defaultBcSys <- c("SP-400", "SP-500", "SP-650", "SP-800", "FP-500", "P-Woolf", "P-KonTiki")
  if (identical(bcSys, "default")) bcSys <- defaultBcSys

  ## Default learning rates (per year) by variant
  # Keep relative ranking unchanged; KonTiki = 0 learning.
  if (is.null(learningRates)) {
    learningRates <- list(low = 0.010, central = 0.005, high = 0)
  }

  ## Read costs data sources
  # Columns expected: biopyrNone, biopyrHeat, biopyrElec, biopyrCHP, biopyrliq
  costRaw <- readSource("PyrolysisCosts", subtype = "REM_BCSys_Costs")

  ## Map keys required for later calculations to property labels in the raw cost source
  propMap <- c(
    inco0_main      = "inco0_main (USD/tDM a)", # per dry matter biomass feedstock
    inco0_low       = "inco0_low (USD/tDM/a)", # per dry matter biomass feedstock
    omf             = "omf (-)",
    omv_total_tDM   = "omv_total_tDM (total USD/tDM)"
  )

  incoMain <- collapseNames(costRaw[, , propMap[["inco0_main"]]])
  incoLow  <- collapseNames(costRaw[, , propMap[["inco0_low"]]])
  omfShare       <- collapseNames(costRaw[, , propMap[["omf"]]])
  omvTotal       <- collapseNames(costRaw[, , propMap[["omv_total_tDM"]]])

  omfShare[is.na(omfShare)] <- 0.04
  omvTotal[is.na(omvTotal)] <- 0

  ## Link plant set-ups to process conditions for yields and biochar energy content
  setupToProcessCond <- c(
    biopyrNone = "SP-500",
    biopyrHeat = "SP-500",
    biopyrElec = "SP-500",
    biopyrCHP  = "SP-500",
    biopyrliq  = "P-Woolf"
  )
  plantSetups <- names(setupToProcessCond)

  ## Get biomass to biochar conversion and biochar attributes
  eff <- collapseNames(calcOutput("BiocharEfficiency", processCond = "default",
                                  outputType = c("dm_yield"), aggregate = FALSE))
  att <- calcOutput("BiocharAttributes", processCond = "default", aggregate = FALSE)

  # dm_yield: tDM biochar / tDM biomass
  # ge: PJ/Mt DM or MJ/kg DM
  # GJ biochar per tDM biomass = dm_yield * ge
  gjOutPerTdm <- collapseNames(eff * att[, , "ge"])

  ## Gate price per plant setup (USD/GJ biochar) in the baseYear (pre-learning)
  gjOutPerTdmPsetup <- new.magpie("GLO", names = plantSetups,
                                  sets = c("region", "year", "plant_setup"))

  for (s in plantSetups) {
    gjOutPerTdmPsetup["GLO", , s] <- gjOutPerTdm["GLO", , setupToProcessCond[[s]]]
  }

  ## Select from data source the initial investment cost according to variant
  incoUseCentral <- incoMain
  incoUseLow <- incoLow
  miss <- is.na(incoUseLow)
  incoUseLow[miss] <- incoMain[miss] * fallbackMultLow
  incoUseHigh <- incoMain

  ## Conversion-only cost per tDM input (USD/tDM input) per variant
  # include initial investment cost, variable operation and maintenance cost,
  # an fixed operation adn maintenance cost
  costCentral <- incoUseCentral * (capRecF + omfShare) + omvTotal
  costLow     <- incoUseLow     * (capRecF + omfShare) + omvTotal
  costHigh    <- (incoUseHigh   * (capRecF + omfShare) + omvTotal) * (1 + markupHigh)

  ## Gate price per plant setup: USD/GJ biochar
  gateLow <- add_dimension(costLow / gjOutPerTdmPsetup, dim = 3.2, add = "variant", nm = "low")
  gateCentral <- add_dimension(costCentral / gjOutPerTdmPsetup, dim = 3.2, add = "variant", nm = "central")
  gateHigh <- add_dimension(costHigh / gjOutPerTdmPsetup, dim = 3.2, add = "variant", nm = "high")

  gatePriceSetupVar <- mbind(gateLow, gateCentral, gateHigh)


  ## Aggregate plant-set ups to bc_tech63 (pyrolysis conditions) per variant
  gatePriceBiocharSysVar <- new.magpie("GLO", names = defaultBcSys,
                                       sets = c("region", "year", "bc_tech63"))
  gatePriceBiocharSysVar <- add_dimension(gatePriceBiocharSysVar, dim = 3.2,
                                          add = "variant", nm  = c("low", "central", "high"))

  ## Direct mapping: P-Woolf from biopyrliq
  gatePriceBiocharSysVar["GLO", , "P-Woolf"] <- gatePriceSetupVar["GLO", , "biopyrliq"]


  ## Aggregate SP-500 from the four plant-setups setups using this pyrolysis condition
  spSetupsAll <- c("biopyrNone", "biopyrHeat", "biopyrElec", "biopyrCHP")
  spSetupsCentral <- c("biopyrElec", "biopyrCHP")

  # For the low variant, the lowest cost across all plant set-ups is used
  gatePriceBiocharSysVar[, , "SP-500"][, , "low"] <-
    min(gatePriceSetupVar[, , spSetupsAll][, , "low"], na.rm = TRUE)

  # For the high variant, the highest cost across all plant set-ups is used
  gatePriceBiocharSysVar[, , "SP-500"][, , "high"] <-
    max(gatePriceSetupVar[, , spSetupsAll][, , "high"], na.rm = TRUE)

  # For the central variant, the mean across all biopyrElec and biopyrCHP is used
  gatePriceBiocharSysVar[, , "SP-500"][, , "central"] <-
    mean(gatePriceSetupVar[, , spSetupsCentral][, , "central"], na.rm = TRUE)


  ## Fill remaining bc_tech63 by scaling from SP-500
  # Keep only systems that are part of bc_tech63
  sysCostScaling <- sysCostScaling[names(sysCostScaling) %in% defaultBcSys]

  for (s in names(sysCostScaling)) {
    # Check if this bc system is still empty
    if (all(!is.finite(as.numeric(gatePriceBiocharSysVar[, , s])))) {
      gatePriceBiocharSysVar[, , s] <- collapseNames(gatePriceBiocharSysVar[, , "SP-500"]) * sysCostScaling[[s]]
    }
  }

  ## Subset requested bc systems
  gatePriceBiocharSysVar <- gatePriceBiocharSysVar[, , intersect(bcSys, defaultBcSys)]



  ## Select on which system level the output should be returned:
  # "bc_sys" (default) returns gate prices aggregated to bc_tech63
  # "biopyr" returns gate prices for individual REMIND plant set-ups

  out <- NULL
  if (outputSys == "biopyr") {
    out  <- gatePriceSetupVar
  } else if (outputSys == "bc_sys") {
    out  <- gatePriceBiocharSysVar
  }


  ## Introduce time dimension and apply learning if selected
  shortYears <- findset("t_all")

  yrsNum <- as.numeric(gsub("^y", "", shortYears))
  baseNum <- as.numeric(gsub("^y", "", baseYear))
  dt <- pmax(0, yrsNum - baseNum)   # 0 for years <= baseYear

  outTime <- new.magpie("GLO", years = shortYears,
                        names = getNames(out),
                        sets  = getSets(out))


  if (is.logical(learningRates)) {
    # If no learning is assumed, costs stay constant over time

    if (!learningRates) {
      for (t in getYears(outTime)) {
        outTime[, t, ] <- setYears(out, t)
      }
    }
  } else {
    # If learning is assumed, costs change over time

    for (v in variant) {

      for (sys in getNames(out, dim = 1)) {
        baseCost <- out[, , v][, , sys]
        lRates <- learningRates[[v]]

        #no learning is assumed for Kontiki
        if (sys == "P-KonTiki") {
          lRates <- 0
        }
        mult <- (1 - lRates)^dt

        for (i in seq_along(shortYears)) {
          outTime[, shortYears[i], sys][, , v] <- baseCost * mult[i]
        }
      }
    }
  }


  return(list(
    x = outTime,
    weight = NULL,
    unit = "USD17MER per GJ biochar",
    description = paste0(
      "Conversion-only levelized biochar gate prices (feedstock excluded) ",
      "for low, central, and high cost assumptions"
    )
  ))
}
