#' @title calcNINWaste
#' @description
#' Calculates the ratio between food supply at household level and food intake for different MAgPIE commodities
#' accounting for food-specific estimates of baseline intake of quantification of NIN diets by the
#' NIN comission, as well as for FAO food waste shares.
#'
#' @param out_type ratio: total food supply to total intake.
#' ratio_detailed_calib: calibrated food-specific estimates.
#' ratio_detailed: food-specific estimates based on FAO food waste shares
#' calib: factor for calibrating estimates based on FAO waste shares to food supply
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link{calcOutput}}, \code{\link{readNIN}},
#' \code{\link{calcNINDiets}}, \code{\link{convertNIN}}
#' @examples
#' \dontrun{
#' calcOutput("NINWaste")
#' }
#' @export

calcNINWaste <- function(out_type = "ratio") {

  fsupply.hist <- calcOutput(type = "FoodSupplyPast", aggregate = FALSE, per_capita = TRUE, product_aggr = FALSE, attributes = "kcal")[, "y2010", ]
  getSets(fsupply.hist)[3] <- "kfo"

  Mag_Intake <- calcOutput("Intake", modelinput = "TRUE", standardize = FALSE, method = "FAO_WHO_UNU1985", aggregate = FALSE)
  Mag_Intake <- collapseNames(Mag_Intake[, "y2010", "SSP2"])

  Mag_NIN_diets <- calcOutput(type = "NINDiets", aggregate = FALSE, attributes = c("kcal", "wm"), calib = TRUE, FAOcountr = FALSE)
  Mag_NIN_diets <- collapseNames(Mag_NIN_diets[, "y2010", "BMK"][, , "2100kcal"][, , "wm", invert = TRUE])

  Intake_calib <- Mag_Intake / dimSums(Mag_NIN_diets, dim = 3)
  Intake_calib[which(!is.finite(Intake_calib))] <- 1

  Mag_NIN_diets <- Mag_NIN_diets * Intake_calib


  #### Calculation of the ratio between food supply at household level and food intake
  # based on FAO estimates on food waste at consumption level and food conversion factors

  FAO_waste_shr <- readSource(type = "FAOLossesWaste", subtype = "Consumption")

  # mapping of FAO waste categories to MAgPIE food commodities
  Mag_kfo <- findset("kfo")

  FAO_wgroups <- c("Oilseeds and pulses", "Cereals", "Roots and tubers", "Fish and seafood", "Oilseeds and pulses", "MNIN", "Eggs", "Milk", "MNIN", "MNIN",
                   "Cereals", "Roots and tubers", "Oilseeds and pulses", "Fruits and vegetables", "Roots and tubers", "Oilseeds and pulses", "Oilseeds and pulses", "Cereals",
                   "Oilseeds and pulses", "Oilseeds and pulses", "Oilseeds and pulses", "Roots and tubers", "Roots and tubers", "Oilseeds and pulses", "Cereals", "Cereals")
  rel_matrix <- cbind(Mag_kfo, FAO_wgroups)

  FAO_waste_shr_detailed <- toolAggregate(FAO_waste_shr, rel = rel_matrix,
                                        dim = 3, from = "FAO_wgroups", to = "Mag_kfo", partrel = FALSE)

  # Conversion factors into edible matter: 0.82 for roots, 0.79 for maize, 0.78 for whNIN, 1 for rice,
  # 0.78 for other grains, 0.77 for fruits and vegetables, 1 for mNIN, 1 for oilseeds and pulses, 1 for milk
  conv_fact <- dimSums(FAO_waste_shr, dim = 1)
  conv_fact[, , ] <- 1
  conv_fact[, , "Cereals"] <- 0.78
  conv_fact[, , "Roots and tubers"] <- 0.82
  conv_fact[, , "Oilseeds and pulses"] <- 1
  conv_fact[, , "Fruits and vegetables"] <- 0.77
  conv_fact[, , "MNIN"] <- 1
  conv_fact[, , "Milk"] <- 1

  conv_fact_detailed <- toolAggregate(conv_fact, rel = rel_matrix,
                                    dim = 3, from = "FAO_wgroups", to = "Mag_kfo", partrel = FALSE)

  conv_fact_detailed[, , "brans"] <- 1
  conv_fact_detailed[, , "maiz"] <- 0.79
  conv_fact_detailed[, , "rice_pro"] <- 1


  # calculation of the ratio of available food at household level to intake (accounting also for losses from food conversion into edible matter)
  overcons_FAO <- 1 / (1 - FAO_waste_shr_detailed) / conv_fact_detailed

  fsupply_estimated <- Mag_NIN_diets * overcons_FAO
  fsupply_calib <-  dimSums(fsupply.hist, dim = 3) / dimSums(fsupply_estimated, dim = 3)
  fsupply_calib[which(!is.finite(fsupply_calib))] <- 1

  overcons_calib <- overcons_FAO * fsupply_calib


  Mag_overcons_fctr <- dimSums(fsupply.hist, dim = 3) / Mag_Intake
  if (min(Mag_overcons_fctr) == 0) {
    temp_overcons <- Mag_overcons_fctr
    temp_overcons[which(Mag_overcons_fctr == 0)] <- NA
    replacement <- as.magpie(apply(temp_overcons, 3, mean, na.rm = TRUE))
    Mag_overcons_fctr[which(Mag_overcons_fctr == 0)] <- replacement
  }

  if (out_type == "ratio") {
    data.out <- Mag_overcons_fctr
    description <- "ratio between total calorie supply and total calorie intake"
  } else if (out_type == "ratio_detailed_calib") {
    data.out <- overcons_calib
    description <- "food-specific ratio between calorie supply and intake consistent with NIN-Lancet baseline diets"
  } else if (out_type == "ratio_detailed") {
    data.out <- overcons_FAO
    description <- "food-specific ratio between calorie supply and intake based on FAO food waste shares and conversion factors"
  } else if (out_type == "calib") {
    data.out <-  fsupply_calib
    description <- "factor for calibrating estimated food supply (based on NIN-Lancet baseline diets and FAO waste shares) to FAO food supply"
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
