#' @title calcForestFireLoss
#'
#' @description Calculate how much loss of forest area happens due to fire disturbances based on FRA 2020 data
#'
#' @return MAgPIE object with FRA 2020 forest fire area loss
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{readFRA2020}}
#' 
#' @examples
#' 
#' \dontrun{ 
#' calcOutput("ForestFireLoss",aggregate=FALSE)
#' }
#' 
#' 

calcForestFireLoss <- function(){
  
    fao <-setNames(readSource("FRA2020",subtype = "forest_fire",convert = TRUE),"overall") # Convert=T returns Mha area lost to fire
    fao <- dimSums(fao,dim=2)/length(getYears(fao)) # Average forest loss due to fires
    curtis <- readSource("ForestLossDrivers")
    
    out <- mbind(fao,curtis)
  
  return(list(x=out,
              weight=NULL,
              unit="Mha",
              description="Area affected by fires in forests",
              isocountries=FALSE))
}
