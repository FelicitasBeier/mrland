#' @title calcEATLancetWaste
#' @description
#' Calculates the ratio between food supply at household level and food intake for different MAgPIE commodities
#' accounting for food-specific estimates of baseline intake of quantification of EAT Lancet diets by the
#' EAT-Lancet comission, as well as for FAO food waste shares.
#'
#' @param out_type ratio: total food supply to total intake.
#' ratio_detailed_calib: calibrated food-specific estimates.
#' ratio_detailed: food-specific estimates based on FAO food waste shares
#' calib: factor for calibrating estimates based on FAO waste shares to food supply
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link{readEATLancet}},
#' \code{\link{calcEATLancetDiets}}, \code{\link{convertEATLancet}}
#' @examples
#' \dontrun{
#' calcOutput("EATLancetWaste")
#' }
calcEATLancetWaste <- function(out_type = "ratio") { # nolint: object_name_linter.

  fsupplyHist <- calcOutput(type = "FoodSupplyPast", aggregate = FALSE, per_capita = TRUE,
                            product_aggr = FALSE, attributes = "kcal")[, "y2010", ]
  getSets(fsupplyHist)[3] <- "kfo"

  magIntake <- calcOutput("Intake", modelinput = TRUE, standardize = FALSE, method = "FAO_WHO_UNU1985",
                          aggregate = FALSE)
  magIntake <- collapseNames(magIntake[, "y2010", "SSP2"])

  magEATDiets <- calcOutput(type = "EATLancetDiets", aggregate = FALSE, attributes = c("kcal", "wm"),
                            calib = TRUE, FAOcountr = FALSE)
  magEATDiets <- collapseNames(magEATDiets[, "y2010", "BMK"][, , "2100kcal"][, , "wm", invert = TRUE])

  intakeCalib <- magIntake / dimSums(magEATDiets, dim = 3)
  intakeCalib[which(!is.finite(intakeCalib))] <- 1

  magEATDiets <- magEATDiets * intakeCalib


  #### Calculation of the ratio between food supply at household level and food intake
  # based on FAO estimates on food waste at consumption level and food conversion factors

  faoWasteShr <- readSource(type = "FAOLossesWaste", subtype = "Consumption")

  # mapping of FAO waste categories to MAgPIE food commodities
  magKfo <- findset("kfo")

  faoWgroups <- c("Oilseeds and pulses", "Cereals", "Roots and tubers", "Fish and seafood",
                  "Oilseeds and pulses", "Meat", "Eggs", "Milk", "Meat", "Meat",
                  "Cereals", "Roots and tubers", "Oilseeds and pulses", "Fruits and vegetables",
                  "Roots and tubers", "Oilseeds and pulses", "Oilseeds and pulses", "Cereals",
                  "Oilseeds and pulses", "Oilseeds and pulses", "Oilseeds and pulses", "Roots and tubers",
                  "Roots and tubers", "Oilseeds and pulses", "Cereals", "Cereals")
  relMatrix <- cbind(magKfo = magKfo, faoWgroups = faoWgroups)

  faoWasteShrDetailed <- toolAggregate(faoWasteShr, rel = relMatrix,
                                       dim = 3, from = "faoWgroups", to = "magKfo", partrel = FALSE)

  # Conversion factors into edible matter: 0.82 for roots, 0.79 for maize, 0.78 for wheat, 1 for rice,
  # 0.78 for other grains, 0.77 for fruits and vegetables, 1 for meat, 1 for oilseeds and pulses, 1 for milk
  convFact <- dimSums(faoWasteShr, dim = 1)
  convFact[, , ] <- 1
  convFact[, , "Cereals"] <- 0.78
  convFact[, , "Roots and tubers"] <- 0.82
  convFact[, , "Oilseeds and pulses"] <- 1
  convFact[, , "Fruits and vegetables"] <- 0.77
  convFact[, , "Meat"] <- 1
  convFact[, , "Milk"] <- 1

  convFactDetailed <- toolAggregate(convFact, rel = relMatrix,
                                    dim = 3, from = "faoWgroups", to = "magKfo", partrel = FALSE)

  convFactDetailed[, , "brans"] <- 1
  convFactDetailed[, , "maiz"] <- 0.79
  convFactDetailed[, , "rice_pro"] <- 1


  # calculation of the ratio of available food at household level to intake
  # (accounting also for losses from food conversion into edible matter)
  overconsFAO <- 1 / (1 - faoWasteShrDetailed) / convFactDetailed

  fsupplyEstimated <- magEATDiets * overconsFAO
  fsupplyCalib <-  dimSums(fsupplyHist, dim = 3) / dimSums(fsupplyEstimated, dim = 3)
  fsupplyCalib[which(!is.finite(fsupplyCalib))] <- 1

  overconsCalib <- overconsFAO * fsupplyCalib


  magOverconsFctr <- dimSums(fsupplyHist, dim = 3) / magIntake
  if (min(magOverconsFctr) == 0) {
    tempOvercons <- magOverconsFctr
    tempOvercons[which(magOverconsFctr == 0)] <- NA
    replacement <- as.magpie(apply(tempOvercons, 3, mean, na.rm = TRUE))
    magOverconsFctr[which(magOverconsFctr == 0)] <- replacement
  }

  if (out_type == "ratio") {
    dataOut <- magOverconsFctr
    description <- "ratio between total calorie supply and total calorie intake"
  } else if (out_type == "ratio_detailed_calib") {
    dataOut <- overconsCalib
    description <- "food-specific ratio between calorie supply and intake consistent with EAT-Lancet baseline diets"
  } else if (out_type == "ratio_detailed") {
    dataOut <- overconsFAO
    description <- paste0("food-specific ratio between calorie supply and intake based on ",
                          "FAO food waste shares and conversion factors")
  } else if (out_type == "calib") {
    dataOut <-  fsupplyCalib
    description <- paste0("factor for calibrating estimated food supply ",
                          "(based on EAT-Lancet baseline diets and FAO waste shares) to FAO food supply")
  } else {
    stop("unknown type")
  }

  dataOut <- setYears(dataOut, NULL)

  #### Define weights and units

  weightPop <- collapseNames(calcOutput("Population", scenario = "SSP2", aggregate = FALSE)[, "y2010", ])
  unit <- "-"

  return(list(x = dataOut,
              weight = weightPop,
              unit = unit,
              description = description))
}
