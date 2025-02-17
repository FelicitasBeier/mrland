#' @title calcFAOLossesWaste
#' @description
#' Calculates the ratio between food supply at household level and food intake for different MAgPIE commodities
#' based on estimated/assumed FAO waste shares for each commodity group (optionally also including food
#' conversion factors into edible matter).
#'
#' @param out_type waste: food-specific ratios based on FAO food waste shares
#' waste_edible:  food-specific ratios based on FAO food waste shares including conversion into edible matter
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link{calcOutput}}, \code{\link{readFAOLossesWaste}},
#' \code{\link{calcEATLancetWaste}}
#' @examples
#' \dontrun{
#' calcOutput("FAOLossesWaste")
#' }
#' @export

calcFAOLossesWaste <- function(out_type = "waste") {

  FAO_waste_shr <- readSource(type = "FAOLossesWaste", subtype = "Consumption")

  # mapping of FAO waste categories to MAgPIE food commodities
  Mag_kfo <- findset("kfo")

  FAO_wgroups <- c("Oilseeds and pulses", "Cereals", "Roots and tubers", "Fish and seafood", "Oilseeds and pulses", "Meat", "Eggs", "Milk", "Meat", "Meat",
                   "Cereals", "Roots and tubers", "Oilseeds and pulses", "Fruits and vegetables", "Roots and tubers", "Oilseeds and pulses", "Oilseeds and pulses", "Cereals",
                   "Oilseeds and pulses", "Oilseeds and pulses", "Oilseeds and pulses", "Roots and tubers", "Roots and tubers", "Oilseeds and pulses", "Cereals", "Cereals")
  rel_matrix <- cbind(Mag_kfo, FAO_wgroups)

  FAO_waste_shr_detailed <- toolAggregate(FAO_waste_shr, rel = rel_matrix,
                                        dim = 3, from = "FAO_wgroups", to = "Mag_kfo", partrel = FALSE)

  # Conversion factors into edible matter: 0.82 for roots, 0.79 for maize, 0.78 for wheat, 1 for rice,
  # 0.78 for other grains, 0.77 for fruits and vegetables, 1 for meat, 1 for oilseeds and pulses, 1 for milk
  conv_fact <- dimSums(FAO_waste_shr, dim = 1)
  conv_fact[, , ] <- 1
  conv_fact[, , "Cereals"] <- 0.78
  conv_fact[, , "Roots and tubers"] <- 0.82
  conv_fact[, , "Oilseeds and pulses"] <- 1
  conv_fact[, , "Fruits and vegetables"] <- 0.77
  conv_fact[, , "Meat"] <- 1
  conv_fact[, , "Milk"] <- 1

  conv_fact_detailed <- toolAggregate(conv_fact, rel = rel_matrix,
                                    dim = 3, from = "FAO_wgroups", to = "Mag_kfo", partrel = FALSE)

  conv_fact_detailed[, , "brans"] <- 1
  conv_fact_detailed[, , "maiz"] <- 0.79
  conv_fact_detailed[, , "rice_pro"] <- 1



  if (out_type == "waste") {
    data.out <- 1 / (1 - FAO_waste_shr_detailed)
    description <- "food-specific ratio between calorie supply and intake based on FAO food waste shares"
  } else if (out_type == "waste_edible") {
    data.out <-  1 / (1 - FAO_waste_shr_detailed) / conv_fact_detailed
    description <- "food-specific ratio between calorie supply and intake based on FAO food waste shares including conversion into edible matter"
  } else {
stop("unknown type")
}

  data.out <- setYears(data.out, NULL)



  #### Define weights and units

  weight.pop <- collapseNames(calcOutput("Population", scenario = "SSP2", aggregate = FALSE)[, "y2010", ])
  unit <- "-"

  return(list(x = data.out,
              weight = weight.pop,
              unit = unit,
              description = description)
  )
}
