#' @title calcEATFruitvegRatio
#' @description
#' Calculates the share of fruits and vegetables in the calorie supply from the others MAgPIE commodity for the past.
#' Information on the calorie supply from fruits and vegetables is relevant in the context of dietary recommendations,
#' e.g. as proposed by the EAT-Lancet Commission on healthy diets from sustainable food systems.
#'
#' @param populationweight datasource of populationweight: FAO can be selected in order to better meet exact values.
#' Normal datasource is PopulationPast
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Isabelle Weindl, Felicitas Beier
#' @seealso \code{\link[madrat]{calcOutput}}, \code{\link{calcEATLancetTargets}},
#' \code{\link[mrfaocore]{calcFAOharmonized}}, \code{\link{calcEATLancetDiets}}
#' @examples
#' \dontrun{
#' calcOutput("EATFruitvegRatio")
#' }
#' @importFrom madrat toolFillWithRegionAvg
#' @importFrom magpiesets findset
#' @export

calcEATFruitvegRatio <- function(populationweight = "PopulationPast") {
  ### FAO Commodity balance
  cbsFAO  <- calcOutput(type = "FAOharmonized", src = "join2010", aggregate = FALSE)[, , "food_supply_kcal"]
  pastYrs <- findset("past_til2020")
  cbsFAO  <- collapseNames(cbsFAO[, pastYrs, ])
  getSets(cbsFAO) <- c("region", "year", "ItemCodeItem")

  ### Population weight
  if (populationweight == "PopulationPast") {
    weight <- collapseNames(calcOutput("PopulationPast", aggregate = FALSE))
  } else if (populationweight == "FAO") {
    weight <- collapseNames(readSource(type = "FAO", subtype = "Pop",
                                       convert = TRUE)[, , "population"]) / 1e+6
  }
  weight <- weight[, pastYrs, ]

  ### Sectoral mapping for FAO items
  relationmatrix <- toolGetMapping("FAOitems_online_2010update.csv", type = "sectoral", where = "mrfaocore")
  relationmatrix <- relationmatrix[, c("post2010_FoodBalanceItem", "k")]
  relationmatrix <- relationmatrix[!duplicated(relationmatrix[, "post2010_FoodBalanceItem"]), ]

  ### Definitions of the food groups "others" and "fruitvegOthers"
  ### (all fruits and vegetables contained in "others")
  others <- relationmatrix[relationmatrix$k %in% "others", "post2010_FoodBalanceItem"]
  fruitvegOthers <- c("2601|Tomatoes and products",
                      "2602|Onions",
                      "2605|Vegetables, other",
                      "2611|Oranges, Mandarines",
                      "2612|Lemons, Limes and products",
                      "2613|Grapefruit and products",
                      "2614|Citrus, Other",
                      "2617|Apples and products",
                      "2618|Pineapples and products",
                      "2619|Dates",
                      "2620|Grapes and products (excl wine)",
                      "2625|Fruits, other")
  ### Definitions of the food groups "cassava" and "bananas"
  ### (bananas and plantains contained in "cassav_sp")
  cassava <- relationmatrix[relationmatrix$k %in% "cassav_sp", "post2010_FoodBalanceItem"]
  bananas <- c("2615|Bananas",
               "2616|Plantains")

  # check if all fruits and vegetables included in "fruitvegOthers" are also contained in "others":
  if (any(setdiff(fruitvegOthers, others)) == TRUE) {
    stop("The above definition of fruits and vegetables (fruitvegOthers)
         is not consistent with the current mapping of FAO items to MAgPIE commodities.")
  }
  # check if all fruits included in "bananas" are also contained in "cassava":
  if (any(setdiff(bananas, cassava)) == TRUE) {
    stop("The above definition of fruits (bananas and plantains)
         is not consistent with the current mapping of FAO items to MAgPIE commodities.")
  }

  ### Aggregation of calories and calculation of ratios:
  # For fruits in others category
  kcalOthers   <- dimSums(cbsFAO[, , others], dim = 3)
  kcalFruitVeg <- dimSums(cbsFAO[, , fruitvegOthers], dim = 3)

  ratioFruitveg2Others <- kcalFruitVeg / kcalOthers
  ratioFruitveg2Others <- toolFillWithRegionAvg(ratioFruitveg2Others, weight = weight, verbose = FALSE)
  ratioFruitveg2Others <- add_dimension(ratioFruitveg2Others, add = "EAT_special", nm = "others")

  # For fruits in cassav_sp category
  kcalCassava <- dimSums(cbsFAO[, , cassava], dim = 3)
  kcalBananas <- dimSums(cbsFAO[, , bananas], dim = 3)

  ratioFruit2Cassava <- kcalBananas / kcalCassava
  # Correct NA's in the past by extending values from years that have values
  if (any(is.na(ratioFruit2Cassava))) {
    noCassavaRegions <- where(is.na(ratioFruit2Cassava))$true$regions
    for (i in noCassavaRegions) {
      yrsNoNA <- where(!is.na(ratioFruit2Cassava[i, , ]))$true$years
      if (!identical(yrsNoNA, character(0))) {
        tmp <- new.magpie(cells_and_regions = i,
                          years = yrsNoNA)
        tmp[, , ] <- ratioFruit2Cassava[i, , ][!is.na(ratioFruit2Cassava[i, , ])]
        ratioFruit2Cassava[i, , ] <- toolHoldConstant(tmp, pastYrs)
      }
    }
  }
  # Fill missing countries with regional average
  ratioFruit2Cassava <- toolFillWithRegionAvg(ratioFruit2Cassava, weight = weight, verbose = FALSE)
  ratioFruit2Cassava <- add_dimension(ratioFruit2Cassava, add = "EAT_special", nm = "cassav_sp")

  out <- mbind(ratioFruitveg2Others, ratioFruit2Cassava)

  min <- 0

  return(list(x = out,
              weight = weight,
              unit = "-",
              description = paste0("share of fruits and vegetables in the others food group ",
                                   "and cassava_sp food group respectively"),
              min = min))
}
