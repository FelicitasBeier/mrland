#' @title calcForestDisturbances
#'
#' @description Calculates which share of forest land is lost due to forest disturbances (including insects, diseases, severe weather events and other causes)
#'
#' @return MAgPIE object with FRA 2020 forest disturbance shares
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{readFRA2020}}
#' 
#' @examples
#' 
#' \dontrun{ 
#' calcOutput("ForestDisturbances",aggregate=FALSE)
#' }
#' 
#' 

calcForestDisturbances <- function(){
  
  disturbance_area   <- setNames(readSource("FRA2020",subtype = "disturbance",convert = TRUE),NULL) # Convert=T returns Mha area lost to fire
  forest_area <- setNames(readSource("FRA2020",subtype = "forest_area",convert = TRUE)[,,"naturallyRegeneratingForest"],NULL) # Convert=T returns Mha area of forests
  
  common_yrs <- intersect(getYears(disturbance_area),getYears(forest_area))
  
  disturbance_share <- disturbance_area[,common_yrs,]/forest_area[,common_yrs,]
  disturbance_share[is.infinite(disturbance_share)] <- 0
  disturbance_share[is.na(disturbance_share)] <- 0
  disturbance_share[disturbance_share>1] <- 1
  
  out <- setYears(disturbance_share[,"y2010",],NULL)
  weight <- setYears(forest_area[,"y2010",],NULL)
  
  return(list(x=out,
              weight=weight,
              unit="1",
              description="Area affected by forest disturbances including insects, diseases, severe weather events and other causes",
              isocountries=FALSE))
}
