#' @title calcForestFireShare
#'
#' @description Calculates which share of forest land is lost due to forest fires
#'
#' @return MAgPIE object with FRA 2020 forest fire shares
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{readFRA2020}}
#' 
#' @examples
#' 
#' \dontrun{ 
#' calcOutput("ForestFireShare",aggregate=FALSE)
#' }
#' 
#' 

calcForestFireShare <- function(){
  
  fire_area   <- setNames(readSource("FRA2020",subtype = "forest_fire",convert = TRUE),NULL) # Convert=T returns Mha area lost to fire
  forest_area <- setNames(readSource("FRA2020",subtype = "forest_area",convert = TRUE)[,,"naturallyRegeneratingForest"],NULL) # Convert=T returns Mha area of forests
  
  common_yrs <- intersect(getYears(fire_area),getYears(forest_area))
  
  fire_share <- fire_area[,common_yrs,]/forest_area[,common_yrs,]
  fire_share[is.infinite(fire_share)] <- 0
  fire_share[is.na(fire_share)] <- 0
  fire_share[fire_share>1] <- 1
  
  out <- setYears(fire_share[,"y2000",],NULL)
  weight <- setYears(forest_area[,"y2000",],NULL)
  
  return(list(x=out,
              weight=weight,
              unit="1",
              description="Area affected by fires in forests",
              isocountries=FALSE))
}
