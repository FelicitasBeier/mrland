#' @title downloadInocencio2007
#'
#' @description Create (no external download) unit cost table based on manual
#'              extraction from Inocencio et al. (2007)
#'
#' @return A list with metadata required by \code{downloadSource()}.
#'
#' @author Felicitas Beier
#'
#' @seealso \code{\link[madrat]{downloadSource}}
#' @examples
#' \dontrun{downloadSource("Inocencio2007)
#' }
#'
#' @export
#'
downloadInocencio2007 <- function() {

  # Data manually extracted from PDF document
  # Table 7: Average unit costs of irrigation projects by region, by purpose of
  # project and by success and failure case (in USD/ha (in 2000 prices))
  # Average over period 1965-1999
  # Success projects only
  unitCost <- data.frame(region = c("GLO", "SSA", "MEA", "SA", "SEA", "EA", "LAM"),
                         new_totalCost       = c(4785, 5726, 8464, 2526, 3861, 4101, 3663),
                         new_hardwareCost    = c(3748, 3556, 7044, 2141, 3146, 3294, 2841),
                         rehab_totalCost     = c(1969, 3488, 3193,  898,  965, 1990, 3730),
                         rehab_hardwareCost  = c(1488, 2303, 2383,  674,  711, 1735, 3004),
                         stringsAsFactors = FALSE)

  # ---- write to file in current working directory (set by downloadSource) ----
  fname <- "unitCost_Inocencio2007.csv"
  utils::write.csv(unitCost, file = fname, row.names = FALSE)

  # ---- metadata for madrat wrapper ----
  meta <- list(url         = "https://www.researchgate.net/publication/42764931_Costs_and_Performance_of_Irrigation_Projects_A_Comparison_of_Sub-Saharan_Africa_and_Other_Developing_Regions",
               title       = "Costs and Performance of Irrigation Projects: A Comparison of Sub-Saharan Africa and Other Developing Regions",
               author      = "Inocencio, A.; Kikuchi, M.; Tonosaki, M.; Maruyama, A.; Merrey, D.; Sally, H.; de Jong, I.",
               license     = "Copyright by IWMI. All rights reserved",
               ISBN        = "978-92-9090-658-2",
               description = paste0("Region-level unit costs for new irrigation development and rehabilitation ",
                                    "based on Table 7 of IWMI Research Report 109 by Inocencio et al. (2007). ",
                                    "Selected only success projects. ",
                                    "Extracted following regions: Sub-Saharan Africa (SSA), Middle East and North Africa (MEA), ",
                                    "South Asia (SA), Southeast Asia (SEA), East Asia (EA), Latin America and the Caribbean (LAM), ",
                                    "and assigned average over all success projects across these regions as rest of the world average (GLO)."),
              unit        = "USD/ha (in 2000 prices)",
              reference   = "Inocencio et al. (2007). Costs and Performance of Irrigation Projects: A Comparison of Sub-Saharan Africa and Other Developing Regions. IWMI Research Report 109.")
  return(meta)
}
