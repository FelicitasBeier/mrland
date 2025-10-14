#' @title readStrefler2021
#' @description Reads in a reporting mif file from REMIND
#' 
#' @param subtype Either "intensive" or "extensive"
#' @return MAgPIE object with regional aggregation of REMIND H12
#' @author Florian Humpenöder
#' @seealso
#' \code{\link[madrat]{readSource}}
#' @examples
#' 
#' \dontrun{ 
#' readSource("Strefler2021",aggregate=FALSE)
#' }
#' @importFrom magclass read.report

readStrefler2021 <- function(subtype) {
  
#https://www.nature.com/articles/s41467-021-22211-2
  
  file_list <- c("CarbonPriceTrajectories.csv")
  
  if(grepl("_",subtype)){
    
    subtype     <- strsplit(subtype, split="_")
    revision    <- as.numeric(unlist(subtype)[2])

  }
  
  x <- NULL
  for(f in file_list) {
    x <- mbind(x,read.report(f,as.list = FALSE))
  }
  
  # remove model name
  x <- collapseNames(x)
  
  #select scenarios and set names
  sel <- c("NPI","exptax","exptax_peak","budget_peak","hotelling2const","GDPgrowth2const","taxLinear_1300")
  names(sel) <- paste("PIK",c("NPI","HOS","HBL","OPT","H2C","GDP","LIN"),sep="_")
  
  x <- x[,,sel]
  getNames(x,dim=1) <- names(sel)
  
  return(x)

}
