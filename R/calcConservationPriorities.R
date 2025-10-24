#' @title calcConservationPriorities
#' @description Function calculates land area in conservation priority areas that was unprotected in 2020 (WDPA).
#'
#' @param consvBaseYear Reference year for land conservation. Chosing "y1750", for instance, means that
#' the reference land use is based on the year 1750 ('pre-industrial') so
#' land use can be restored to the pre-industrial state in conservation priority areas.
#' Any year available in the LUH3 data set can be chosen. Historic land use in the LUH3 data
#' is based on the HYDE data base.
#' The choice "y2020" provides a special case, in which reference land use is based on the 2020 ESA CCI LULC map,
#' derived at a spatial resolution of 300 x 300 Meter.
#' @param cells (deprecated) number of cells of landmask (only "lpjcell" for 67420 cells)
#' @param nclasses Options are either "seven" or "nine".
#' \itemize{
#' \item "seven" separates primary and secondary forest and includes "crop", "past", "forestry",
#' "primforest", "secdforest", "urban" and "other"
#' \item "nine" adds the separation of pasture and rangelands, as well as a differentiation of
#' primary and secondary non-forest vegetation and therefore returns "crop", "past", "range",
#' "forestry", "primforest", "secdforest", "urban", "primother" and "secdother"
#' }
#'
#' @return magpie object in cellular resolution with different protection options in conservation priority areas
#' @author Patrick v. Jeetze
#'
#' @examples
#' \dontrun{
#' calcOutput("ConservationPriorities", aggregate = FALSE)
#' }
#'
#' @importFrom magpiesets findset addLocation
#' @importFrom magclass collapseDim
#' @importFrom mstools toolCoord2Isocell
#'

calcConservationPriorities <- function(consvBaseYear = "y1750", cells = "lpjcell", nclasses = "seven") {
  # ===============================
  # Get conservation templates
  # ===============================

  # ------------------------------
  # Key Biodiversity Areas (KBAs)
  # ------------------------------
  # KBAs have global value for conservation, due to their
  # outstanding ecological integrity, globally important ecosystems
  # or significant populations of animals, fungi and plants.
  # https://www.iucn.org/resources/conservation-tool/key-biodiversity-areas

  kba <- calcOutput("KeyBiodiversityAreas",
    maginput = TRUE, unprotected = TRUE,
    nclasses = nclasses, aggregate = FALSE
  )

  # ---------------------------
  # Global Safety Net (GSN)
  # ---------------------------
  # Areas of the terrestrial realm where increased
  # conservation action is needed to protect biodiversity
  # and store carbon (Dinerstein et al. 2020, Science).

  gsn <- calcOutput("GlobalSafetyNet",
    maginput = TRUE, nclasses = nclasses,
    aggregate = FALSE
  )
  # Unprotected Key Biodiversity Areas have been masked at
  # high resolution to remove overlaps and can be now
  # added up
  gsn <- gsn + kba
  gsn <- collapseDim(gsn, dim = 3.3)

  # ----------------------------------
  # Critical Connectivity Areas (CCA)
  # ----------------------------------
  # Brennan et al. (2022), Science:
  # "Areas where the flow of animal movement
  # is concentrated are places with the potential to
  # disproportionately reduce connectivity if further
  # restricted or destroyed [...]"

  # Critical Connectivity Areas relate to th upper
  # 90th percentile of values of mapped
  # mammal movement probabilities

  cca <- calcOutput("CriticalConnectivityAreas",
    maginput = TRUE, nclasses = nclasses,
    mask = "KBA", aggregate = FALSE
  )

  # Unprotected Key Biodiversity Areas have been masked at
  # high resolution to remove overlaps and can be now
  # combined with KBAs to create a balanced conservation
  # template
  cca <- cca + kba
  cca <- collapseDim(cca, dim = 3.3)

  # ------------------------------
  # Irrecoverable carbon
  # ------------------------------
  # Noon et al. (2022), Nature Sustainability:
  # "Irrecoverable carbon reserves are [...] manageable, vulnerable to
  # disturbance and could not be recovered by 2050 if lost today.
  # While irrecoverabilitycan be considered over any timeframe, we
  # selected 30 years as the most policy-relevant scenario to align with
  # the Paris Agreement goal to reach net-zero emissions by mid-century."

  irrC <- calcOutput("IrrecoverableCarbonLand",
    maginput = TRUE, nclasses = nclasses,
    aggregate = FALSE
  )

  # -------------------------------------------------
  # Biodiversity Hotspots & Intact Forest Landscapes
  # -------------------------------------------------
  # Biodiversity Hotspots have at least 1500 vascular plants as endemics
  # and must have 30 % or less of its original natural vegetation.
  # https://www.conservation.org/priorities/biodiversity-hotspots

  # An intact forest landscape (IFL) is a seamless mosaic of forest and
  # naturally treeless ecosystems with no remotely detected signs of
  # human activity and a minimum area of 500 km2 (Potapov et al. 2017,
  # Science Advances)

  bhifl <- calcOutput("BHIFL",
    nclasses = nclasses,
    aggregate = FALSE
  )

  # -------------------------------------------------
  # Half Earth (PBL)
  # -------------------------------------------------
  # An ecoregion-based approach to protecting half of
  # the terrestrial surface. 50 % of land are conserved in
  # each to the 867 ecoregions.

  pblHalfEarth <- calcOutput("HalfEarth",
    nclasses = nclasses,
    aggregate = FALSE
  )

  getNames(pblHalfEarth) <- paste0("PBL_", getNames(pblHalfEarth))

  # ----------------------------
  # 30 % of land surface
  # ----------------------------
  # A scenario for the 30 by 30 target that is based on
  # Key Biodiversity Areas, Distinct Species Assemblages (DSA)
  # from the Global Safety Net and Critical Connectivity Areas (CCAs).

  gsn30 <- calcOutput("GlobalSafetyNet",
    maginput = TRUE, nclasses = nclasses,
    aggregate = FALSE
  )
  # Use GSN Distinct Species Assemblages (DSA) for the 30 % target
  gsn30 <- gsn30[, , "GSN_DSA"]

  # To remove overlaps CCAs have already been masked
  # by KBAs and GSN DSAs at high resolution
  cca30 <- calcOutput("CriticalConnectivityAreas",
    maginput = TRUE, nclasses = nclasses,
    mask = "KBA_GSN", aggregate = FALSE
  )

  # combine KBA, GSN & CCA data sets
  thirty <- mbind(kba, gsn30, cca30)
  thirty <- dimSums(thirty, dim = 3.1)

  # Unprotected KBAs, DSAs and CCAs make up 16.13 %
  # of the land surface. Together with the WDPA
  # protected area in 2020 they make up ~30% of the
  # global land surface.

  getNames(thirty) <- paste("30by30", getNames(thirty), sep = ".")

  # --------------------------------------------------
  # Half Earth (Global Safety Net)
  # --------------------------------------------------
  # The sum of all conservation clusters of the GSN
  # pertains to 50.4 % of the global terrestrial realm
  # together with currently protected areas.

  gsnHalfEarth <- calcOutput("GlobalSafetyNet",
    maginput = TRUE, nclasses = nclasses,
    aggregate = FALSE
  )

  gsnHalfEarth <- mbind(kba, gsnHalfEarth)
  gsnHalfEarth <- dimSums(gsnHalfEarth, dim = 3.1)
  getNames(gsnHalfEarth) <- paste("GSN_HalfEarth", getNames(gsnHalfEarth), sep = ".")

  # ======================================
  # Prepare output
  # ======================================
  # Combine all templates for the output
  consvPrio <- mbind(thirty, kba, gsn, bhifl, irrC, cca, gsnHalfEarth, pblHalfEarth)

  # ----------------------------
  # Correct data
  # ----------------------------
  # Due to small mismatches conservation priority land
  # in the additive options can be larger than total
  # land in a grid cell. This is corrected in the following.

  luh3 <- calcOutput("LUH3", cellular = TRUE, yrs = consvBaseYear,
                     landuseTypes = "LUH3", aggregate = FALSE)[, consvBaseYear, ]
  getYears(luh3) <- NULL

  # get total land area
  totLand <- dimSums(luh3, dim = 3)

  # urban land
  urbanLand <- calcOutput("UrbanLandFuture",
    subtype = "LUH3", aggregate = FALSE,
    timestep = "5year"
  )

  # baseline protected area (WDPA)
  wdpaLand <- calcOutput("ProtectedAreaBaseline",
    aggregate = FALSE, nclasses = "seven",
    magpie_input = TRUE,
  )

  # make sure that future conservation priority land is not greater
  # than total land area minus urban area and currently protected areas
  landNoUrban <- setYears(totLand, "y2020") -
    setCells(urbanLand[, "y2020", "SSP2"], getCells(totLand)) -
    setCells(dimSums(wdpaLand[, "y2020", ], dim = 3), getCells(totLand))
  landNoUrban[landNoUrban < 0] <- 0

  getYears(landNoUrban) <- getYears(consvPrio)
  # compute mismatch factor
  landMismatch <- setNames(landNoUrban, NULL) / dimSums(consvPrio, dim = 3.2)
  landMismatch <- toolConditionalReplace(landMismatch, c(">1", "is.na()"), 1)
  # correct conservation priority data
  consvPrio <- consvPrio * landMismatch

  # -------------------------------
  # Define conservation base year
  # -------------------------------
  # Historic land use is derived from the LUH3 data.
  # Based on historic land use, shares of the different
  # land types are extracted and applied proportionally
  # in conservation priority areas.

  if (consvBaseYear != "y2020") {
    # Reclassify LUH classes to MAgPIE classes
    if (nclasses == "seven") {
      consvBaseLand <- mbind(
        setNames(dimSums(luh3[, , c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")], dim = 3), "crop"),
        setNames(dimSums(luh3[, , c("pastr", "range")], dim = 3), "past"),
        setNames(luh3[, , c("primf", "secdf")], c("primforest", "secdforest")),
        setNames(new.magpie(getCells(luh3), fill = 0), "forestry"),
        luh3[, , c("urban")],
        setNames(dimSums(luh3[, , c("primn", "secdn")], dim = 3), "other")
      )
    } else if (nclasses == "nine") {
      consvBaseLand <- mbind(
        setNames(dimSums(luh3[, , c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")], dim = 3), "crop"),
        setNames(luh3[, , "pastr"], "past"),
        setNames(luh3[, , "range"], "range"),
        setNames(luh3[, , c("primf", "secdf")], c("primforest", "secdforest")),
        setNames(new.magpie(getCells(luh3), fill = 0), "forestry"),
        luh3[, , c("urban")],
        setNames(luh3[, , "primn"], "primother"),
        setNames(luh3[, , "secdn"], "secdother")
      )
    }

    # calculate share of land types in conservation base year
    consvBaseLandShr <- consvBaseLand / totLand
    consvBaseLandShr <- toolConditionalReplace(consvBaseLandShr, conditions = "is.na()", 0)

    # Multiply total conservation priority land with
    # share of land types in conservation base year
    consvPrio <- dimSums(consvPrio, dim = 3.2) * consvBaseLandShr
  }

  return(list(x = consvPrio,
              weight = NULL,
              unit = "Mha",
              description = paste0("Land conservation priority targets in each land type. ",
                                   "Land use in conservation priority areas is based on the ",
                                   "reference year ", consvBaseYear,
                                   ifelse(gsub("y", "", consvBaseYear) <= 1800, " (pre-industrial)", "")),
              isocountries = FALSE))
}
