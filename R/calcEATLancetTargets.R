#' @title calcEATLancetTargets
#' @description
#' Calculates minimum and maximum targets for healthy food intake according to reference
#' recommendations proposed by the EAT-Lancet Commission on healthy diets from sustainable
#' food systems, specified for different MAgPIE commodities.
#'
#' @param  attributes Attributes of food commodities (available: kcal/d and g/d)
#'
#' @return List of magpie objects with results on country level, weight on country level,
#'         unit and description.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link{calcOutput}}, \code{\link{readEATLancet}},
#' \code{\link{calcEATLancetDiets}}
#' @examples
#' \dontrun{
#' calcOutput("EATLancetTargets")
#' }
#' @export

calcEATLancetTargets <- function(attributes = "kcal/d") {

  eatTargets <- readSource(type = "EATLancet", subtype = "recommend")
  eatTargets[which(is.na(eatTargets))] <- 0

  # mapping of EAT Lancet food categories to MAgPIE food commodities
  eatFood    <- c("legumes", "fish",
                  "poultry", "eggs", "milk",
                  "sugar", "oil_palm", "oil_veg",
                  "nuts_seeds",
                  "vegetables", "fruits",
                  "fruits_starch",   "roots",   "fg_redmeat")
  eatTFood15 <- c("t_legumes", "t_fish",
                  "t_livst_chick", "t_livst_egg", "t_livst_milk",
                  "t_sugar", "t_oils", "t_oils",
                  "t_nutseeds",
                  "t_fruitveg", "t_fruitveg",
                  "t_fruitstarch", "t_roots", "t_redmeat")

  relMatrixFood <- cbind(eatFood, eatTFood15)

  eat <- toolAggregate(eatTargets, rel = relMatrixFood,
                       dim = 3.1, partrel = TRUE,
                       from = "eatFood", to = "eatTFood15")


  #### Selection of relevant attributes
  out <- collapseNames(eat[, , attributes])

  #### Define weights and units
  weight <- collapseNames(calcOutput("Population", scenario = "SSP2", aggregate = FALSE)[, "y2010", ])
  unit   <- "kcal or grams wm per capita per day"

  return(list(x = out,
              weight = weight,
              unit = unit,
              description = paste0("Minimum and maximum targets for healthy food intake ",
                                   "according to the EAT-Lancet Commission"))
  )
}
