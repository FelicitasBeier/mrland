#' @title calcPlantationContribution
#' @description
#' Calculates the interpolated contribution share of plantations to roundwood demand
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' \dontrun{
#' calcOutput("PlantationContribution")
#' }
#' @importFrom magclass getNames<- as.magpie time_interpolate
#' @importFrom stringr str_extract_all
#' @importFrom madrat toolGetMapping
#' @export

calcPlantationContribution <- function() {
  ## Read Share from source
  out <- mbind(readSource("TimberShare", subtype = "abare", convert = TRUE),
               readSource("TimberShare", subtype = "brown", convert = TRUE))

  out <- setYears(out, "y1995")

  ## Create Scenarios
  scen <- c("constant", "h5s5l5", "h5s2l2", "h5s2l1", "h5s1l1", "h5s1l05", "h2s1l05")

  ## Years
  year <- paste0("y", seq(2000, 2250, 5))

  hist  <- paste0("y", seq(1995, 2020, 5))
  short <- paste0("y", seq(2025, 2050, 5))
  long  <- paste0("y", seq(2055, 2250, 5))

  out <- time_interpolate(dataset = out, interpolated_year = year, integrate_interpolated_years = TRUE,
                          extrapolation_type = "constant")
  out <- add_dimension(x = out, dim = 3.2, nm = scen)

  outScen <- out[, , "constant", invert = TRUE]

  mods <- str_extract_all(scen, "\\d+")
  names(mods) <- scen

  for (i in getNames(outScen, dim = "new")) {
    for (j in 2:length(year)) {
      present <- getYears(outScen)[j]
      past    <- getYears(outScen)[j - 1]

      scenPattern <- mods[[i]]
      multiplier <- c((1 + as.numeric(mods[[i]]) / 100))

      if (length(grep(pattern = "0", x = scenPattern, value = TRUE)) > 0) {
        pos <- match(grep(pattern = "0", x = scenPattern, value = TRUE), scenPattern)
        multiplier[pos] <- (1 + as.numeric(grep(pattern = "0", x = scenPattern, value = TRUE)) / 1000)
      }

      if (present %in% hist) outScen[, present, i] <- setYears(outScen[, past, i], NULL) * multiplier[1]
      if (present %in% short) outScen[, present, i] <- setYears(outScen[, past, i], NULL) * multiplier[2]
      if (present %in% long) outScen[, present, i] <- setYears(outScen[, past, i], NULL) * multiplier[3]
    }
  }
  out[, , getNames(outScen)] <- outScen
  out[, length(year) + 1, ] <- out[, length(year), ]
  out[out > 1] <- 1
  out <- round(out, 3)

  ## Fix JPN and REF values - Both have 0 plantations Growing stock so set
  ## Plantation contribution to overall demand in JPN = 0 (or very tiny value)
  ## Plantation contribution to overall demand in REF = 0.01

  ## Find standard mapping - which countries belong to REF and JPN in standard mapping -
  ## we will modify ISO codes here so that this works with all mappings
  h12mapping <- toolGetMapping(type = "regional", name = "regionmappingH12.csv", where = "madrat")
  jpn <- h12mapping[h12mapping$RegionCode == "JPN", ]$CountryCode
  ref <- h12mapping[h12mapping$RegionCode == "REF", ]$CountryCode
  eur <- h12mapping[h12mapping$RegionCode == "EUR", ]$CountryCode
  out[jpn, , ] <- 0.00001
  out[ref, , ] <- 0.01
  out[eur, , ] <- out[eur, , ] * 3
  ## WARNING : (DO NOT CHANGE 0.01 value in REF to any value lower than this as
  ## this will result in wrong plantation establishment in REF and would need
  ## adjustment in establishment calibration factors - currently 0.05 for REF)

  out[out > 1] <- 1

  ## Weight
  weight <- collapseNames(calcOutput("TimberDemand", aggregate = FALSE)[, "y1995", "production"])[, , c("Roundwood")]

  return(list(x = out,
              weight = weight,
              min = 0,
              unit = "percent",
              description = "Calculates the share of roundwood production coming from timber plantations"))
}
