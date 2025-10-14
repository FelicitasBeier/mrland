#' @title calcEndUseTimber
#' @description 
#' Calculates the demand of timber from historical FAO data (including intermediate products).
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Abhijeet Mishra
#' @seealso
#' \code{\link[mrfaocore]{calcFAOmassbalance_pre}}
#' @examples
#' 
#' \dontrun{ 
#' calcOutput("EndUseTimber")
#' }
#' @importFrom magpiesets findset 
#' @export

calcEndUseTimber <- function(){

  fao_timber <- collapseNames(calcOutput("TimberDemand",aggregate = F)[,,"domestic_supply"])
  
  fao_timber <- fao_timber[,c(paste0("y",seq(from=1965,to=2015,by=5))),]
  
  fao_timber[,"y2015","Wood fuel"] <- fao_timber[,"y2010","Wood fuel"]
  
  getNames(fao_timber) <- gsub(x = getNames(fao_timber),pattern = " ",replacement = "_")
  getNames(fao_timber) <- gsub(x = getNames(fao_timber),pattern = "-",replacement = "_")
  getNames(fao_timber) <- tolower(getNames(fao_timber))
  
  out <- fao_timber
  
  return(list(x=out,weight=NULL,
                unit="mio m3/yr",
                description="historical timber domestic demand")
         )
}
