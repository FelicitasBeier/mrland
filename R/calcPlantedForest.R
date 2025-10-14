#' @title calcPlantedForest
#' @description
#' Calculates the share of plantations in planted forest
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' \dontrun{
#' calcOutput("PlantedForest")
#' }
#' @importFrom magclass getNames<- as.magpie time_interpolate
#' @export

calcPlantedForest <- function() {
  ## Read land area frpom source
  a <- readSource("FRA2020", "forest_area")
  plantedShare <- setNames(round(a[, , "plantationForest"], 3) / round(a[, , "plantedForest"], 3), NULL)
  plantedShare[is.na(plantedShare)] <- 0
  out <- setYears(plantedShare[, "y2000", ], NULL)

  ## Change EUR values - See Forestry GMD paper review from Pekka Lauri
  magIsoReg <- toolGetMapping(type = "regional", name = "regionmappingH12.csv", where = "madrat")
  regEur <- magIsoReg$CountryCode[magIsoReg$RegionCode == "EUR"]
  out[regEur, , ] <- out[regEur, , ] * 3
  out[regEur, , ][out[regEur, , ] > 1] <- 1

  ## Weight
  weight <- setYears(setNames(round(a[, "y2000", "plantedForest"], 3), NULL), NULL)

  return(list(x = out,
              weight = weight,
              min = 0,
              unit = "share",
              description = "Calculates the share of plantation forest in planted forest"))

}
