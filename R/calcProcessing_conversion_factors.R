#' @title calcProcessing_conversion_factors
#'
#' @description Calculates global conversion factors
#' from the FAOmassbalance_pre in order to feed to MAgPIE for
#' conversion of primary to secondary products
#'
#' @return magpie object of conversion factors
#' @author Benjamin Bodirsky, David M Chen
#' @examples
#' \dontrun{
#' a <- calcOutput("Processing_conversion_factors",
#'                 aggregate = FALSE)
#' }
#' @importFrom magclass setNames getNames

calcProcessing_conversion_factors <- function() { # nolint: object_name_linter.
  massbalance <- calcOutput("FAOmassbalance", aggregate = FALSE)

  kpr <- findset("kpr")
  ksd <- findset("ksd")
  kprocessing <- findset("processing20")
  mbReduced <- dimSums(massbalance[, , "dm"], dim = c(1, 3.3))
  kmb <- paste("X", kpr, sep = "")

  productionEstimated <- dimSums(mbReduced[, , "production_estimated"][, , ksd], dim = 3.2)
  convmatrix <- add_dimension(x = productionEstimated, dim = 3.2, add = "kpr", nm = kmb)
  convmatrix <- add_dimension(x = convmatrix, dim = 3.1, add = "processing", nm = kprocessing)
  convmatrix[, , ] <- 0

  mbReduced2 <- mbReduced[, , kpr]

  tmp <- dimSums(mbReduced2[, , c("alcohol1", "alcohol2", "alcohol3", "alcohol4")], dim = 3.2) /
    dimSums(mbReduced2[, , "fermentation"], dim = 3.2)
  convmatrix[, , "alcohol"][, , "fermentation"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "brewers_grain1"], dim = 3.2) /
    dimSums(mbReduced2[, , "fermentation"], dim = 3.2)
  convmatrix[, , "distillers_grain"][, , "fermentation"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "brans1"], dim = 3.2) /
    dimSums(mbReduced2[, , "milling"], dim = 3.2)
  convmatrix[, , "brans"][, , "milling"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "branoil1"], dim = 3.2) /
    dimSums(mbReduced2[, , "milling"], dim = 3.2)
  convmatrix[, , "oils"][, , "milling"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "oilcakes1"], dim = 3.2) /
    dimSums(mbReduced2[, , "milling"], dim = 3.2)
  convmatrix[, , "oilcakes"][, , "milling"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "distillers_grain1"], dim = 3.2) /
    dimSums(mbReduced2[, , "distilling"], dim = 3.2)
  convmatrix[, , "distillers_grain"][, , "distilling"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "ethanol1"], dim = 3.2) /
    dimSums(mbReduced2[, , "distilling"], dim = 3.2)
  convmatrix[, , "ethanol"][, , "distilling"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "molasses1"], dim = 3.2) /
    dimSums(mbReduced2[, , "refining"], dim = 3.2)
  convmatrix[, , "molasses"][, , "refining"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , c("sugar1", "sugar2", "sugar3")], dim = 3.2) /
    dimSums(mbReduced2[, , "refining"], dim = 3.2)
  convmatrix[, , "sugar"][, , "refining"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , "oilcakes1"], dim = 3.2) /
    dimSums(mbReduced2[, , "extracting"], dim = 3.2)
  convmatrix[, , "oilcakes"][, , "extracting"] <- setNames(tmp, paste0("X", getNames(tmp)))

  tmp <- dimSums(mbReduced2[, , c("oil1", "oil2")], dim = 3.2) /
    dimSums(mbReduced2[, , "extracting"], dim = 3.2)
  convmatrix[, , "oils"][, , "extracting"] <- setNames(tmp, paste0("X", getNames(tmp)))

  # add conversion attributes of Single cell Protein (SCP)
  # based on Table S3 in Pikaar et al 2018, which provides conversion factors
  # as "ton substrate DM / ton microbial protein (MP)"
  # Miscanthus (begr): 5.5 t DM begr / t DM MP
  # Sugar Cane (sugr_cane): 4.3 t DM sugar_cane / t DM MP
  # maize (foddr): 5.6 t DM foodr / t DM MP
  # We need to convert these factors into "Conversion factors of primary products into secondary products"
  # How much MP (secondary product) do get per ton DM substrate (primary product)?
  # For each t DM sugr_cane we get 1/4.3=0.2326 t DM MP.
  convmatrix[, , "breeding"][, , "scp"][, , "Xbegr"] <- 1 / 5.5 # 0.1818
  convmatrix[, , "breeding"][, , "scp"][, , "Xsugr_cane"] <- 1 / 4.3 # 0.2326
  convmatrix[, , "breeding"][, , "scp"][, , "Xfoddr"] <- 1 / 5.6 # 0.1786

  tmp <- collapseNames(dimSums(massbalance[, , "production"][, , "dm"][, , "fibres"], dim = 1) /
                         dimSums(massbalance[, , "production"][, , "dm"][, , "cottn_pro"], dim = 1))
  convmatrix[, , "ginning"][, , "fibres"][, getYears(tmp), "Xcottn_pro"] <- tmp
  convmatrix[is.nan(convmatrix)] <- 0
  convmatrix[is.na(convmatrix)] <- 0
  convmatrix[is.infinite(convmatrix)] <- 0
  # add conversion attributes of bioenergy crops into ethanol

  convmatrix[, , "extracting"][, , "ethanol"][, , "Xbetr"] <- 0.36
  convmatrix[, , "extracting"][, , "ethanol"][, , "Xbegr"] <- 0.36

  # test that the conversion factors are congruent with the
  # production_estimated in FAOmassbalance_pre
  kprocessingM <- setdiff(kprocessing, c("breeding", "ginning"))
  mbReduced <- dimSums(massbalance[, , "dm"], dim = 3.3)

  a <- mbReduced[, , kprocessingM][, , kpr]
  getNames(a, dim = 1) <- paste0("X", getNames(a, dim = 1))
  b <- dimSums(a * convmatrix[, , kprocessingM], dim = c(3.1, 3.2))
  c <- dimSums(mbReduced[, , ksd][, , "production_estimated"], dim = 3.2)
  if (any(round(b - c, 4) != 0)) {
    stop("conversion factors are incoherent with production_estimated column of massbalance!")
  }

  getNames(convmatrix, dim = 3) <- substring(getNames(convmatrix, dim = 3), 2)

  convmatrix <- toolHoldConstantBeyondEnd(convmatrix)

  return(list(x = convmatrix,
              weight = NULL,
              unit = "t DM / t DM",
              description = paste0("Conversion factors of primary products into secondary products. ",
                                   "primary product x conversion factor = secondary product")))
}
