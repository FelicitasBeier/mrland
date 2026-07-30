#' @title calcClimateRegionsIPCC
#' @description calculates IPCC Climate Regions (IPCC2006 ch.4.3) based on t, ppt,
#' pet from LPJml. elevation dimension not included for tropical montane class
#'
#' @param landusetypes all or only one
#' @param cellular FALSE for country level, TRUE for cells
#' @param convert fills missing countries for country level aggregation with warm
#' temperate moist (mostly small island nations)
#' @param yearly FALSE for normal magpie 5 year time spans, TRUE for yearly
#'
#' @return Country or cellular magpie object with fraction of each climate region by country or cell
#' @author David Chen
#'
#' @examples
#' \dontrun{
#' calcOutput("ClimateRegionsIPCC")
#' }
#' @importFrom SPEI thornthwaite
#' @importFrom magpiesets findset
#' @importFrom mrlandcore toolLPJmLDefault

calcClimateRegionsIPCC <- function(landusetypes = "all", cellular = FALSE,
                                   yearly = FALSE, convert = TRUE) {
  # extract default arguments for LPJmL
  cfg <- toolLPJmLDefault()

  # PET based on thornwaite function for historical baseline
  pet <- calcOutput("LPJmLTransform", lpjmlversion = cfg$defaultLPJmLVersion,
                    climatetype = cfg$baselineHist,
                    subtype     = "pnv:pet", subdata = NULL,
                    stage       = "smoothed:cut",
                    monthly     = TRUE,
                    aggregate   = FALSE)

  # LPJmL climate inputs: precipitation and temperature
  precip <- calcOutput("LPJmLClimateInput", climatetype = cfg$baselineHist,
                       variable = "precipitation:monthlySum", stage = "smoothed",
                       lpjmlVersion = cfg$defaultLPJmLVersion, aggregate = FALSE)

  temp <- calcOutput("LPJmLClimateInput", climatetype = cfg$baselineHist,
                     variable = "temperature:monthlyMean", stage = "smoothed",
                     lpjmlVersion = cfg$defaultLPJmLVersion, aggregate = FALSE)

  # yearly total precip and pet
  precip <- dimSums(precip, dim = c(3.1))
  pet    <- dimSums(pet, dim = c(3.1))
  # yearly MAP:PET ratio
  ratio <- precip / (pet + 0.000000001)
  # binary where ratio>1
  ratio <- ratio > 1
  # Mean Annual Temperature
  temp <- dimSums(temp, dim = 3.1) / 12
  # if any month temp is above 10, binary TRUE
  t <- temp > 10
  tMonthly10 <- dimSums(t, dim = 3) > 0

  getNames(temp)   <- "temp"
  getNames(ratio)  <- "ratio"
  getNames(precip) <- "precip"
  getNames(tMonthly10) <- "t_monthly10"
  clm <- mbind(temp, precip, ratio, tMonthly10)
  # climate zones, set to zero, note, no tropical montane climate (exists in IPCC classification)
  clm <- add_columns(clm, dim = 3.1, addnm = c("tropical_wet", "tropical_moist", "tropical_dry",
                                               "warm_temp_moist", "warm_temp_dry",
                                               "cool_temp_moist", "cool_temp_dry",
                                               "boreal_moist", "boreal_dry",
                                               "polar_moist", "polar_dry"))
  clm[, , c("tropical_wet", "tropical_moist", "tropical_dry",
            "warm_temp_moist", "warm_temp_dry",
            "cool_temp_moist", "cool_temp_dry",
            "boreal_moist", "boreal_dry",
            "polar_moist", "polar_dry")] <- 0
  # IPCC decision tree
  clm[, , "tropical_wet"][which(clm[, , "temp"] > 18 & clm[, , "precip"] > 2000)] <- 1
  clm[, , "tropical_moist"][which(clm[, , "temp"] > 18 & clm[, , "precip"] <= 2000 & clm[, , "precip"] > 1000)] <- 1
  clm[, , "tropical_dry"][which(clm[, , "temp"] > 18 & clm[, , "precip"] <= 1000)] <- 1
  clm[, , "warm_temp_moist"][which(clm[, , "temp"] > 10 & clm[, , "temp"] <= 18 & clm[, , "ratio"] == 1)] <- 1
  clm[, , "warm_temp_dry"][which(clm[, , "temp"] > 10 & clm[, , "temp"] <= 18 & clm[, , "ratio"] == 0)] <- 1
  clm[, , "cool_temp_moist"][which(clm[, , "temp"] > 0 & clm[, , "temp"] <= 10 & clm[, , "ratio"] == 1)] <- 1
  clm[, , "cool_temp_dry"][which(clm[, , "temp"] > 0 & clm[, , "temp"] <= 10 & clm[, , "ratio"] == 0)] <- 1
  clm[, , "boreal_moist"][which(clm[, , "temp"] < 0 & clm[, , "t_monthly10"] == 1 & clm[, , "ratio"] == 1)] <- 1
  clm[, , "boreal_dry"][which(clm[, , "temp"] < 0 & clm[, , "t_monthly10"] == 1 & clm[, , "ratio"] == 0)] <- 1
  clm[, , "polar_moist"][which(clm[, , "temp"] < 0 & clm[, , "t_monthly10"] == 0 & clm[, , "ratio"] == 1)] <- 1
  clm[, , "polar_dry"][which(clm[, , "temp"] < 0 & clm[, , "t_monthly10"] == 0 & clm[, , "ratio"] == 0)] <- 1

  climateRegions <- clm[, , c("tropical_wet", "tropical_moist", "tropical_dry",
                              "warm_temp_moist", "warm_temp_dry",
                              "cool_temp_moist", "cool_temp_dry",
                              "boreal_moist", "boreal_dry",
                              "polar_moist", "polar_dry")]

  if (cellular == TRUE) {
    out <- climateRegions
  } else if (cellular == FALSE) {

    if (landusetypes == "all") {
      landuse <- calcOutput("LanduseInitialisation", cells = "lpjcell", cellular = TRUE, aggregate = FALSE)
    } else if (landusetypes != "all") {
      landuse <- calcOutput("LanduseInitialisation", cells = "lpjcell",
                            cellular = TRUE, aggregate = FALSE)[, , landusetypes]
    }
    # cell to country aggregation
    landuse <- dimOrder(collapseDim(addLocation(landuse), dim = c("cell")),  perm = c(2, 3, 1), dim = 1)
    names(dimnames(landuse))[1] <- "x.y.iso"
    landuse <- time_interpolate(landuse, interpolated_year = getYears(ratio), extrapolation_type = "constant")
    landuseSum <- dimSums(landuse, dim = 3)
    climateRegions <- toolAggregate(x = climateRegions, weight = landuseSum, to = "iso")
    out <- climateRegions

    if (convert == TRUE) {
      out <- toolCountryFill(out, fill = 0)
      missing <- setdiff(getItems(out, dim = 1.1), getItems(climateRegions, dim = 1.1))
      out[missing, , "warm_temp_moist"] <- 1
    }
  }

  weight <- calcOutput("LanduseInitialisation", cellular = FALSE, aggregate = FALSE)
  weight <- time_interpolate(weight, interpolated_year = getYears(out), extrapolation_type = "constant")
  if (landusetypes != "all") {
    weight <- weight[, , landusetypes]
  }

  return(list(x = out,
              weight = weight,
              unit = "NULL",
              min = 0,
              max = 1,
              description = "Proportion of IPCC Climate Region",
              isocountries = !cellular))
}
