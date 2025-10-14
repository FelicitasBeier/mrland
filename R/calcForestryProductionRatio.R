#' @title calcForestryProductionRatio
#' @description 
#' Calculates the management factor(s) needed to upscale the yield of forest plantations as compared to natural vegetation based on FAO data.
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' 
#' \dontrun{ 
#' calcOutput("ForestryProductionRatio")
#' }
#' @importFrom magclass getNames<- as.magpie
#' @export

calcForestryProductionRatio <- function(){
  
  ## Define a mapping
  
  ## Reading in the area data from FAO
  prod_ratio <- readSource("ForestryProductionRatio")

  out <- prod_ratio
  
  weight <- out
  weight[weight>=0] <- 1
  
  return(list(x=out,
              weight=weight,
              min=0,
              unit="factor",
              description="Calculates share of timber production coming from forestry."))
  
}
