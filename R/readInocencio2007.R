#' @title readInocencio2007
#'
#' @description reads in average unit cost of irrigation projects per region and
#'              by purpose of project from Inocencio et al. (2007) IWMI Research Report 109
#'
#' @return magpie object of unit costs of irrigation projects
#'
#' @author Felicitas Beier
#'
#' @seealso \code{\link[madrat]{readSource}}
#' @examples
#' \dontrun{ a <- readSource("Inocencio2005)
#' }
#'
readInocencio2007 <- function() {

  # Read in data by Inocencio et al. (2007), table 7:
  # Average unit costs of irrigation projects by region, by purpose of
  # project and by success and failure case (in USD/ha (in 2000 prices))
  # Average over period 1965-1999
  # Success projects only
  unitCost <- read.csv("unitCost_Inocencio2007.csv")

  # Transform to magpie object
  unitCost <- as.magpie(unitCost, spatial = 1, temporal = 0, datacol = 2)

  # Average over period 1965-1999 assigned to MAgPIE initialization year
  getItems(unitCost, dim = 2) <- "y1995"

return(list(x        = unitCost,
            unit     = "USD per hectare (in 2000 prices)",
            metadata = paste0("Inocencio et al. (2007). Costs and Performance of ",
                              "Irrigation Projects: A Comparison of Sub-Saharan Africa ",
                              "and Other Developing Regions. IWMI Research Report 109.")))
}
