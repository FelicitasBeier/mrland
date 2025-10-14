#' @title readWBirrigation
#' 
#' @description reads in World bank irrigation data: WBirrigation data .csv file as magclass object
#' from Jones, William I. 1995. "World Bank and Irrigation." Washington, D.C.: World Bank. 
#' Bonsch et al. (2015) "Environmental Flow Provision: Implications for Agricultural Water and Land-Use at the Global Scale": Table A1 - Investment costs for expanding irrigation infrastructure in US$ per hectare.
#' Based on: World Bank Irrigation Investment Cost Data. 
#' William I. Jones(1991) "The World Bank and Irrigation" (World Bank Operations Evaluation Study)
#' 
#' @return magpie object of the WBirrigation data
#' 
#' @author Lavinia Baumstark
#' 
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' 
#' \dontrun{ a <- readSource(type="WBirrigation")
#' }
#' 
readWBirrigation <- function() {   
   
   irr     <- read.csv("irrigation.csv", sep=";", row.names=1)
   regions <- c("East and South Asia"="ESA","East Asia"="EAS", "South Asia"="SAS", "India"="IND", "Europe"="ERP", "Middle East"="MET", "Africa" ="AFR", "North Africa"="NAF", "Sub-Saharan Africa"="SAF", "Latin America and Caribbean"="LAC")
   irr     <- data.frame(ad_unit_cost=irr$ad_unit_cost, reg=as.character(regions), stringsAsFactors=F)
   irr     <- as.magpie(irr)
   getYears(irr) <- "y1995"
   
   return(irr)
}  
