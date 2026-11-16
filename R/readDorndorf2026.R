#' @title readDorndorf2026
#' @description Reads REMIND reporting mif files from the biochar model runs of
#' Dorndorf et al (2026), stored per scenario family in subfolders.
#'
#' @param subtype Available subtypes: "Baseline"
#' @return Magpie object with biochar production per plant setup. Sets are
#' "bc_plant63" (plus "biopyrelec", which calcBiocharProd adds to "biopyrchp")
#' and "scenario".
#' @author Kristine Karstens
#' @seealso
#' \code{\link[madrat]{readSource}}, \code{\link{readStrefler2021}}
#' @examples
#'
#' \dontrun{
#' readSource("Dorndorf2026", subtype = "Baseline", aggregate = FALSE)
#' }
#' @importFrom magclass read.report

readDorndorf2026 <- function(subtype = "Baseline") {

  files <- c(Baseline = "Baseline/REMIND_generic_PB650-BCpess-CTS01-noFuel-CHP17-EJ100.mif")

  file <- toolSubtypeSelect(subtype, files)

  x <- read.report(file, as.list = FALSE)

  # dim 3 after read.report is scenario.model.variable, dots inside names are
  # escaped to "p" by read.report, so dots can be used as dimension separator.
  full <- getItems(x, dim = 3)
  var  <- sub("^[^.]*\\.[^.]*\\.", "", full)
  scen <- sub("\\..*$", "", full)

  # Mapping between MAgPIE names for pyrolysis plant set-ups and the REMIND
  # reporting keys used in these runs (old key style), where the reported
  # variable name is: SE|Biochar|<key> (EJ/yr). "KonTiki" is not used (values
  # are numerically zero). "biopyrelec" is mapped to "biopyrchp" in
  # calcBiocharProd.
  techMap <- c(biopyrchp  = "SE|Biochar|biopyrCHP (EJ/yr)",
               biopyrhe   = "SE|Biochar|biopyrHeat (EJ/yr)",
               biopyrliq  = "SE|Biochar|biopyrFuel (EJ/yr)",
               biopyronly = "SE|Biochar|biopyrOnly (EJ/yr)",
               biopyrelec = "SE|Biochar|biopyrElec (EJ/yr)")

  if (!all(techMap %in% var)) {
    stop("Missing biochar variables in ", file, ": ",
         paste(setdiff(techMap, var), collapse = ", "))
  }

  sel <- var %in% techMap
  x <- x[, , full[sel]]

  # build dim 3 as plant.scenario, mirroring the REMIND branch in calcBiocharProd
  x <- do.call(mbind, lapply(techMap, function(variable) {
    keep <- var[sel] == variable
    tmp <- x[, , which(keep)]
    getNames(tmp) <- unique(scen[sel][keep])
    add_dimension(tmp, dim = 3.1, nm = names(techMap)[techMap == variable])
  }))

  return(x)
}
