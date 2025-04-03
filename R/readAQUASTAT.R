#' @title readAQUASTAT
#' @description Read in data based on AQUASTAT database
#'              (https://www.fao.org/aquastat/statistics/query/index.html)
#'
#' @param subtype
#' \itemize{
#' \item \code{ConsAgri}:      4454|Conservation agriculture area (1000 ha)
#'                             4454_conservation_agriculture_area_in_1000_ha.csv
#' \item \code{ConsAgriShare}: 4455|Commoditiy Balance LivestockConservation
#'                                  agriculture area as % of arable land area (%)
#'                             4455_conservation_agriculture_area_as_share_of_
#'                                 arable_land_areas.csv)
#' \item \code{rf2irRatio}:    Ratio between rainfed and irrigated yields (%)
#'                             Ratio_between_rainfed_and_irrigated_yields.csv
#' }
#'
#' @return magpie objects with results on contury level
#' @author Kristine Karstens
#' @examples
#' \dontrun{
#' readSource("AQUASTAT", subtype = "ConsAgri", convert = TRUE)
#' }
#'
#' @importFrom utils read.csv

readAQUASTAT <- function(subtype = "ConsAgri") {

  files <-
    c(ConsAgri      = "4454_conservation_agriculture_area_in_1000_ha.csv",
      ConsAgriShare = "4455_conservation_agriculture_area_as_share_of_arable_land_areas.csv",
      rf2irRatio    = "4557_Ratio_between_rainfed_and_irrigated_yields.csv")

  file <- toolSubtypeSelect(subtype, files)

  aquastat   <- read.csv(file, nrows = 600, blank.lines.skip = TRUE,
                         stringsAsFactors = FALSE)

  aquastat$X <- substr(aquastat$X, 1, 3)

  if (subtype %in% c("rf2irRatio")) aquastat <- aquastat[1:399, ]

  aquastat   <- aquastat[aquastat$X != "", ]

  colyears  <- grep("X[0-9]{4}\\.[0-9]{4}", names(aquastat), value = TRUE)
  years     <- paste0("y", c(substr(head(colyears, 1), 2, 5):
                               substr(tail(colyears, 1), 7, 10)))

  mag <- array(NA, dim = c(length(aquastat$X), length(years), 1),
               dimnames = list(aquastat$X, years, subtype))

  for (iso3 in seq_along(aquastat$X)) {

    valueCols <-
      setdiff(which(!is.na(aquastat[iso3, ]) &
                      aquastat[iso3, ] != "" &
                      !grepl("[a-zA-Z]", aquastat[iso3, ])), c(1, 2))

    for (i in valueCols[valueCols %% 2 != 0]) {
      mag[iso3, paste0("y", aquastat[iso3, i]), ] <-
        as.numeric(aquastat[iso3, i + 1])
    }
  }

  if (subtype == "ConsAgri") mag <- mag / 1000 # unit transform 1000 ha <=> 1 Mha

  return(as.magpie(mag))
}
