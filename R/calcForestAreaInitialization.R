#' @title calcForestAreaInitialization
#' @description
#' Calculates the management factor(s) needed to upscale the yield of forest plantations as compared to natural
#' vegetation based on FAO data.
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' \dontrun{
#' calcOutput("ForestAreaInitialization")
#' }
#' @importFrom magclass getNames<- as.magpie
#' @export

calcForestAreaInitialization <- function() {

  ## Reading in the area data from FAO
  forestryArea <- readSource("FRA2015Doc", "plantation_forest")[, "y1995", ]

  ## Reading NatVeg area from FAO
  natvegArea <- setNames(dimSums(mbind(readSource("FRA2015Doc", "secondary_forest")[, "y1995", ],
                                       readSource("FRA2015Doc", "primary_forest")[, "y1995", ])), "natveg_forest")

  zeroCheck <- where(forestryArea == 0 & natvegArea == 0)$true$regions
  if (length(zeroCheck) != 0) {
    message("Countries without forest area detected. Will be manually corrected ---")
    zeroArea <- where(forestryArea == 0 & natvegArea == 0)$true$regions
    forestryArea[zeroArea, , ] <- median(forestryArea[forestryArea != 0])
    natvegArea[zeroArea, , ] <- median(natvegArea[natvegArea != 0]) * 0.25
  }

  zeroCheck <- where(forestryArea == 0 & natvegArea == 0)$true$regions
  if (length(zeroCheck == 0)) {
    message("Countries without forest area filled with median values.")
  }

  out <- setYears(mbind(setNames(forestryArea, "forestry"), setNames(natvegArea, "natveg")), NULL)

  return(list(x = out,
              weight = NULL,
              min = 0,
              unit = "factor",
              description = "Calculates forestry area for preparing plantation management factors in MAgPIE."))

}
