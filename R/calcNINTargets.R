#' @title calcNINTargets
#' @description
#' Calculates minimum and maximum targets for healthy food intake according to reference
#' recommendations proposed by the NIN on healthy diets from sustainable
#' food systems, specified for different MAgPIE commodities.
#'
#' @param  attributes Attributes of food commodities (available: kcal/d and g/d)
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link{calcOutput}}, \code{\link{readNIN}},
#' \code{\link{calcNINDiets}}
#' @examples
#' \dontrun{
#' calcOutput("NINTargets")
#' }
#' @export

calcNINTargets <- function(attributes = "kcal/d") {

  NIN_rtargets <- readSource(type = "NIN", subtype = "recommend")
  NIN_rtargets[which(is.na(NIN_rtargets))] <- 0


  #### Mapping to MAgPIE

  # mapping of NIN  food categories to MAgPIE food commodities
  NIN_food <-    c("legumes", "soybeans", "fish", "poultry", "eggs", "milk", "sugar", "oil_palm", "oil_veg")
  NIN_tfood15 <- c("t_puls_pro", "t_soybean", "t_fish", "t_livst_chick", "t_livst_egg", "t_livst_milk", "t_sugar", "t_oils", "t_oils")

  rel_matrix_food <- cbind(NIN_food, NIN_tfood15)

  NIN_rtargets_food <- toolAggregate(NIN_rtargets, rel = rel_matrix_food,
                                     dim = 3.1, from = "NIN_food", to = "NIN_tfood15", partrel = TRUE)

  # mapping of NIN  food groups to MAgPIE diet target groups
  NIN_group <-    c("nuts_seeds", "vegetables", "fruits", "roots", "fg_redmNIN")
  NIN_tgroup15 <-    c("t_nutseeds", "t_fruitveg", "t_fruitveg", "t_roots", "t_redmNIN")

  rel_matrix_group <- cbind(NIN_group, NIN_tgroup15)

  NIN_rtargets_group <- toolAggregate(NIN_rtargets, rel = rel_matrix_group,
                                      dim = 3.1, from = "NIN_group", to = "NIN_tgroup15", partrel = TRUE)

  NIN <- mbind(NIN_rtargets_food, NIN_rtargets_group)


  #### Selection of relevant attributes
  NIN.out <- collapseNames(NIN[, , attributes])


  #### Define weights and units
  weight.pop <- collapseNames(calcOutput("Population", scenario = "SSP2", aggregate = FALSE)[, "y2010", ])
  unit <- "kcal or grams wm per capita per day"

  return(list(x = NIN.out,
              weight = weight.pop,
              unit = unit,
              description = "Minimum and maximum targets for healthy food intake according to the NIN")
  )
}
