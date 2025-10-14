#' @title calcForestProductionInitialization
#' @description
#' Calculates the management factor(s) needed to upscale the yield of forest plantations as compared to
#' natural vegetation based on FAO data.
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' \dontrun{
#' calcOutput("ForestProductionInitialization")
#' }
#' @importFrom magclass getNames<- as.magpie
#' @export

calcForestProductionInitialization <- function() { # nolint

  ## Calculate overall timber production
  timbProd <- collapseNames(calcOutput("TimberDemand", aggregate = FALSE)[, "y1995", "production"])[, , c("Roundwood")]

  ## Splitting production -- Calculating production share
  prodPlantShare <- setYears(collapseNames(readSource("TimberShare", subtype = "abare")), NULL)

  ## Create dummy magpie object with production from plantations and natveg
  splitProdn <- setNames(timbProd, "forestry")
  splitProdn <- add_columns(splitProdn, addnm = "natveg", dim = 3.1)

  ## Filling the values
  splitProdn[, , "forestry"] <- setNames(timbProd, NULL) * prodPlantShare
  splitProdn[, , "natveg"] <- setNames(timbProd, NULL) - setNames(splitProdn[, , "forestry"], NULL)

  ## Check if a country has 0 area but still reporting production
  fai <- calcOutput("ForestAreaInitialization", aggregate = FALSE) ## Calling the area data

  ## Check 0 area and positive production in plantations
  zeroCheck <- where(fai[, , "forestry"] == 0 & splitProdn[, , "forestry"] > 0)$true$regions
  if (length(zeroCheck != 0)) {
    message("Countries without plantation area but positive production from plantation detected ---")
    fai[zeroCheck, , "forestry"] <- 0.1
  }

  ## Check 0 area and positive production in natveg
  zeroCheck <- where(fai[, , "natveg"] == 0 & splitProdn[, , "natveg"] > 0)$true$regions
  if (length(zeroCheck != 0)) {
    message("Countries without natveg area but positive production from natveg detected ---")
    fai[zeroCheck, , "natveg"] <- 2
  }

  out <- setYears(mbind(setNames(splitProdn[, , "forestry"], "forestry_prod"),
                        setNames(fai[, , "forestry"], "forestry_area"),
                        setNames(splitProdn[, , "natveg"], "natveg_prod"),
                        setNames(fai[, , "natveg"], "natveg_area")), NULL)

  return(list(x = out,
              weight = NULL,
              min = 0,
              unit = "factor",
              description = "Calculates forestry area for preparing plantation management factors in MAgPIE."))

}
