#' @title calcEATLancetDiets
#' @description
#' Calculates daily per capita intake for MAgPIE food commodities that are consistent with diet scenarios
#' developed by the EAT-Lancet Commission on healthy diets from sustainable food systems.
#' The unit is kcal/day per capita or wm/day per capita.
#' Mapping of intake from EAT Lancet to MAgPIE food commodities is done indivudually for the different available units.
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#'
#' @param attributes attributes of different food commodities (available: kcal and wm).
#' @param calib if TRUE, total daily per capita intake for MAgPIE food commodities is calibrated to EAT Lancet total intake.
#' @param FAOcountr if TRUE, estimates for countries not covered in FAOSTAT are set to Zero.
#'
#' @author Isabelle Weindl
#' @seealso \code{\link{calcOutput}}, \code{\link{readEATLancet}},
#' \code{\link{convertEATLancet}}
#' @examples
#' \dontrun{
#' calcOutput("EATLancetDiets")
#' }
#' @export
#' @importFrom magclass getSets<- new.magpie add_dimension dimSums where
#' @importFrom madrat toolNAreplace

calcEATLancetDiets <- function(attributes = c("wm", "kcal"), calib = TRUE, FAOcountr = FALSE) {
  # read data for EAT Lancet diet (intake)
  EATdiets <- readSource(type = "EATLancet", subtype = "cons_data")
  getNames(EATdiets, dim = 2) <- c("wm", "kcal")
  getSets(EATdiets)[4] <- "unit"

  # read data on food supply based on FAOSTAT, aggregated to MAgPIE commodities
  # food supply includes householde waste: food supply = intake + waste
  kfo <- findset("kfo")
  fsupplyHist <- calcOutput(type = "FoodSupplyPast", aggregate = FALSE,
                            per_capita = TRUE, product_aggr = FALSE,
                            attributes = c("wm", "kcal"))
  getSets(fsupplyHist)[3:4] <- c("kfo", "unit")

  # define new diet object with MAgPIE food products
  years <- dimnames(EATdiets)[[2]]
  iso   <- dimnames(EATdiets)[[1]]
  magEATdiets <- new.magpie(cells_and_regions = iso,
                            years = years,
                            names = getNames(EATdiets, dim = 1),
                            sets = c("region", "year", "kcal_scn"))

  magEATdiets <- add_dimension(magEATdiets, add = "unit", dim = 3.2)
  magEATdiets <- add_columns(magEATdiets, addnm = getNames(EATdiets, dim = 2), dim = 3.2)
  magEATdiets <- magEATdiets[, , "dummy", pmatch = TRUE, invert = TRUE]

  magEATdiets <- add_dimension(magEATdiets, add = "diet_scn", dim = 3.3)
  magEATdiets <- add_columns(magEATdiets, addnm = getNames(EATdiets, dim = 3), dim = 3.3)
  magEATdiets <- magEATdiets[, , "dummy", pmatch = TRUE, invert = TRUE]

  magEATdiets <- add_dimension(magEATdiets, add = "kfo", dim = 3.4)
  magEATdiets <- add_columns(magEATdiets, addnm = kfo, dim = 3.4)
  magEATdiets <- magEATdiets[, , "dummy", pmatch = TRUE, invert = TRUE]


  ##### Mapping of EAT Lancet food commodities to MAgPIE commodities:


  # 0. No information needed from EAT Lancet diets since food demand is zero:
  magEATdiets[, , "scp"] <- 0


  # 1. Direct 1:1 mapping possible
  dir_kfo        <- c("fish", "livst_chick", "livst_egg", "livst_milk", "livst_pig",
                      "maiz", "puls_pro", "rice_pro", "soybean", "sugar")
  dir_food_group <- c("fish", "poultry", "eggs", "milk", "pork",
                      "maize", "legumes", "rice", "soybeans", "sugar")
  dir_rel_matrix <- cbind(dir_kfo, dir_food_group)

  magEATdiets[, , dir_kfo] <- toolAggregate(EATdiets[, , dir_food_group], rel = dir_rel_matrix,
                                            from = "dir_food_group", to = "dir_kfo",
                                            dim = 3.4, partrel = FALSE)


  # 2. MAgPIE commodities can be represented by several EAT Lancet commodities
  magEATdiets[, , "livst_rum"] <- EATdiets[, , "beef"] + EATdiets[, , "lamb"]
  magEATdiets[, , "oils"]      <- EATdiets[, , "oil_palm"] + EATdiets[, , "oil_veg"]


  # 3. EAT Lancet commodities are an aggregate of several commodities in MAgPIE

  # roots
  roots <- c("cassav_sp", "potato")
  for (i in roots) {
    root.shr <- setYears(fsupplyHist[, "y2010", i] / dimSums(fsupplyHist[, "y2010", roots],
                                                             dim = 3.1),
                         NULL)
    if (any(!is.finite(root.shr))) {
      replacement <- as.magpie(apply(root.shr, 3, mean, na.rm = TRUE))
      root.shr <- toolNAreplace(root.shr, replaceby = replacement)$x
    }
    magEATdiets[, , i] <- EATdiets[, , "roots"] * root.shr
  }


  # 4. Definition mismatch between MAgPIE and EAT Lancet commodities
  # (not possible to express commodity_dataset_x as sum(commodities_dataset_y))

  ### temperate and tropical cereals
  crls_Mag <- c("tece", "trce")
  crls_EAT <- c("wheat", "othr_grains")
  # estimate waste share for cereals based on 2010 food demand (FAO) and 2010 BMK intake from EAT Lancet
  waste_shr_crls <- 1 - setYears(dimSums(EATdiets[, , "BMK"][, , "2100kcal"][, "y2010", crls_EAT],
                                         dim = 3.4) / dimSums(fsupplyHist[, "y2010", crls_Mag],
                                                              dim = 3.1),
                                 NULL)
  if (any(!is.finite(waste_shr_crls))) {
    temp_waste_shr <- waste_shr_crls
    temp_waste_shr[which(!is.finite(temp_waste_shr))] <- 0
    replacement <- as.magpie(apply(temp_waste_shr, 3, mean, na.rm = TRUE))
    waste_shr_crls <- collapseNames(toolNAreplace(waste_shr_crls, replaceby = replacement)$x)
  }

  mult_factor_othgr <- collapseNames(EATdiets[, , "othr_grains"] / setYears(EATdiets[, "y2010", "othr_grains"][, , "BMK"], NULL))
  if (any(!is.finite(mult_factor_othgr))) {
    temp_mult_factor <- mult_factor_othgr
    temp_mult_factor[which(!is.finite(temp_mult_factor))] <- 1
    replacement <- dimSums(mult_factor_othgr, dim = 1)
    for (t in years) {
      replacement[, t, ] <- as.magpie(apply(temp_mult_factor[, t, ], 3, mean, na.rm = TRUE))
    }
    mult_factor_othgr <- collapseNames(toolNAreplace(mult_factor_othgr, replaceby = replacement)$x)
  }

  balance.post.crls <- collapseNames(EATdiets[, , "othr_grains"] - (
                                                                    setYears(collapseNames(fsupplyHist[, "y2010", "trce"]), NULL) * (1 - waste_shr_crls[, , ]) * mult_factor_othgr))

  balance.post.crls[which(balance.post.crls < 0)] <- 0

  magEATdiets[, , "trce"] <- EATdiets[, , "othr_grains"] - balance.post.crls
  magEATdiets[, , "tece"] <- EATdiets[, , "wheat"] + balance.post.crls


  ### vegetables, fruits and nuts/seeds

  # the EAT Lancet categories "nuts_seeds","vegetables" and "fruits" are mapped
  # to the MAgPIE categories "groundnut","rapeseed","sunflower" and "others"

  vfns_Mag <- c("others", "groundnut", "rapeseed", "sunflower")
  ns_Mag   <- c("groundnut", "rapeseed", "sunflower")
  vfns_EAT <- c("nuts_seeds", "vegetables", "fruits")

  # step 1: estimate "magEATdiets[,,"others"]" by using an intermediate mapping
  #        between EAT Lancet (supercategory "vegfruits" and original "nuts_seeds") and MAgPIE ("others" and supercategory "nutsoilseeds")

  tmp_EAT_vegfruits <- EATdiets[, , "vegetables"] + EATdiets[, , "fruits"]
  tmp_EAT_nuts_seeds <- EATdiets[, , "nuts_seeds"]

  # For the split of EAT Lancet "nuts_seeds" into MAgPIE "others" and "nutsoilseed", we use information on food availability of
  # MAgPIE "others" and "nutsoilseed" and transform it into intake equivalents by using waste assumptions (from FAO) combined with
  # conversion factors into edible matter that were both applied in Springmann et al. (2018) to calculate baseline food consumption
  # from IMPACT food availability data, as explained in the Supplement to:
  # Springmann M, Wiebe K, Mason-D’Croz D, Sulser T, Rayner M, Scarborough P. Health and nutritional aspects of sustainable diet strategies and their
  # association with environmental impacts: a global modelling analysis with country-level detail. Lancet Planet Health 2018

  # Conversion factors into edible matter: 0.77 for fruits and vegetables, 1 for oilseeds and pulses
  conv_fact_vf <- 0.77
  FAO_waste_shr <- readSource(type = "FAOLossesWaste", subtype = "Consumption")

  # calibration of waste (vegetables, fruits, all nuts and seeds), as estimated following Springmann et al. 2019, with Mag_waste = food availability - intake
  # based on information on 2010 food availability in MAgPIE (derived from FAOSTAT) and 2010 BMK intake from EAT Lancet

  Mag_waste_vfns <- collapseNames(dimSums(fsupplyHist[, "y2010", vfns_Mag], dim = 3.1) - dimSums(EATdiets[, , "BMK"][, , "2100kcal"][, "y2010", vfns_EAT], dim = 3.4))

  waste_vfns_estimated <- collapseNames(FAO_waste_shr[, "y2010", "Oilseeds and pulses"] * dimSums(fsupplyHist[, "y2010", ns_Mag], dim = 3.1)
                                        + FAO_waste_shr[, "y2010", "Fruits and vegetables"] * fsupplyHist[, "y2010", "others"]
                                        + (1 - FAO_waste_shr[, "y2010", "Fruits and vegetables"]) * (1 - conv_fact_vf) * fsupplyHist[, "y2010", "others"])

  calib_fact_waste <- Mag_waste_vfns / waste_vfns_estimated

  if (any(!is.finite(calib_fact_waste))) {
    temp_calib <- calib_fact_waste
    temp_calib[which(!is.finite(temp_calib))] <- 1
    replacement <- as.magpie(apply(temp_calib, 3, mean, na.rm = TRUE))
    calib_fact_waste <- collapseNames(toolNAreplace(calib_fact_waste, replaceby = replacement)$x)
  }

  waste_shr_ns <- collapseNames(setYears(calib_fact_waste * FAO_waste_shr[, "y2010", "Oilseeds and pulses"], NULL))

  mult_factor_ns <- collapseNames(tmp_EAT_nuts_seeds / setYears(tmp_EAT_nuts_seeds[, "y2010", "BMK"], NULL))
  if (any(!is.finite(mult_factor_ns))) {
    temp_mult_factor <- mult_factor_ns
    temp_mult_factor[which(!is.finite(temp_mult_factor))] <- 1
    replacement <- dimSums(mult_factor_ns, dim = 1)
    for (t in years) {
      replacement[, t, ] <- as.magpie(apply(temp_mult_factor[, t, ], 3, mean, na.rm = TRUE))
    }
    mult_factor_ns <- collapseNames(toolNAreplace(mult_factor_ns, replaceby = replacement)$x)
  }

  balance.post.vfns <- collapseNames(tmp_EAT_nuts_seeds - (
                                                           setYears(dimSums(fsupplyHist[, "y2010", ns_Mag], dim = 3.1), NULL) * (1 - waste_shr_ns[, , ]) * mult_factor_ns))

  balance.post.vfns[which(balance.post.vfns < 0)] <- 0

  tmp_Mag_nutsoilseeds      <- tmp_EAT_nuts_seeds - balance.post.vfns
  magEATdiets[, , "others"] <- tmp_EAT_vegfruits + balance.post.vfns


  # step 2: now, the intermediate MagPIE supercategory "nutsoilseeds" can be distributed to its components ("groundnut","rapeseed","sunflower"),
  #        similar to the case of roots above

  ns_Mag <- c("groundnut", "rapeseed", "sunflower")

  for (i in ns_Mag) {
    nutsoilseeds.shr <- setYears(fsupplyHist[, "y2010", i] / dimSums(fsupplyHist[, "y2010", ns_Mag], dim = 3.1), NULL)
    if (any(!is.finite(nutsoilseeds.shr))) {
      replacement <- as.magpie(apply(nutsoilseeds.shr, 3, mean, na.rm = TRUE))
      nutsoilseeds.shr <- toolNAreplace(nutsoilseeds.shr, replaceby = replacement)$x
    }
    magEATdiets[, , i] <- tmp_Mag_nutsoilseeds * nutsoilseeds.shr
  }


  # 5. No quantitative information about future contribution of certain foods to diets in EAT Lancet dataset

  ### sugar crop products besides sugar:
  # for the BMK scenario, estimate intake of sugar crop products based on food supply data, but only until the value
  # from the EAT Lancet "othrcrp" category is reached
  magSugrcrp <- c("molasses", "sugr_beet", "sugr_cane")
  magEATdiets[, , magSugrcrp] <- 0

  factorSC <- collapseNames(EATdiets[, , "BMK"][, , "othrcrp"]
                            / setYears(EATdiets[, "y2010", "BMK"][, , "2100kcal"][, , "othrcrp"], NULL))
  if (any(!is.finite(factorSC))) {
    temp_mult_factor <- factorSC
    temp_mult_factor[which(!is.finite(temp_mult_factor))] <- 1
    replacement <- dimSums(factorSC, dim = 1)
    for (t in years) {
      replacement[, t, ] <- as.magpie(apply(temp_mult_factor[, t, ], 3, mean, na.rm = TRUE))
    }
    factorSC <- collapseNames(toolNAreplace(factorSC, replaceby = replacement)$x)
  }

  for (t in years) {
    magEATdiets[, t, magSugrcrp][, , "BMK"] <- setYears(fsupplyHist[, "y2010", magSugrcrp] * (1 - FAO_waste_shr[, "y2010", "Fruits and vegetables"])
                                                        * conv_fact_vf, t) * factorSC[, t, ]

    sugr_gap_reg <- where(dimSums(magEATdiets[, t, magSugrcrp], dim = 3.4) > EATdiets[, t, "BMK"][, , "othrcrp"])$true$regions
    sugr_gap_cal <- collapseNames(EATdiets[sugr_gap_reg, t, "BMK"][, , "2100kcal"][, , "othrcrp"] / dimSums(magEATdiets[sugr_gap_reg, t, magSugrcrp][, , "BMK"][, , "2100kcal"], dim = 3.4))
    sugr_gap_cal[which(sugr_gap_cal > 1)] <- 1

    magEATdiets[sugr_gap_reg, t, magSugrcrp][, , "BMK"] <- collapseNames(magEATdiets[sugr_gap_reg, t, magSugrcrp][, , "BMK"] * sugr_gap_cal[sugr_gap_reg, t, ])
  }

  # estimate the residual amount of "othrcrp" without the amount that is mapped to sugar crops
  res_EAT_othrcrp <- EATdiets[, , "othrcrp"] - dimSums(magEATdiets[, , magSugrcrp], dim = 3.4)
  # mean: 4.3 kcal; min: 0 kcal; max: 231.6 kcal


  ### alcohol:
  # for the BMK scenario, estimate intake of alcohol based on food supply data, but only until the value
  # from the residual EAT Lancet category "othrcrp" plus the EAT Lancet category "othrcal" is reached

  magEATdiets[, , "alcohol"] <- 0

  tmp_EAT_othrcal_pos <- EATdiets[, , "othrcal"]
  tmp_EAT_othrcal_pos[which(tmp_EAT_othrcal_pos < 0)] <- 0
  tmp_EAT_othrcal_neg <- EATdiets[, , "othrcal"]
  tmp_EAT_othrcal_neg[which(tmp_EAT_othrcal_neg > 0)] <- 0

  balance_post_al <- collapseNames(tmp_EAT_othrcal_pos[, , ] + res_EAT_othrcrp)

  mult_factor_al <- collapseNames(balance_post_al[, , "BMK"]
                                  / setYears(balance_post_al[, "y2010", "BMK"][, , "2100kcal"], NULL))
  if (any(!is.finite(mult_factor_al))) {
    temp_mult_factor <- mult_factor_al
    temp_mult_factor[which(!is.finite(temp_mult_factor))] <- 1
    replacement <- dimSums(mult_factor_al, dim = 1)
    for (t in years) {
      replacement[, t, ] <- as.magpie(apply(temp_mult_factor[, t, ], 3, mean, na.rm = TRUE))
    }
    mult_factor_al <- collapseNames(toolNAreplace(mult_factor_al, replaceby = replacement)$x)
  }

  for (t in years) {
    magEATdiets[, t, "alcohol"][, , "BMK"] <- setYears(fsupplyHist[, "y2010", "alcohol"], t) * 0.90 * factorSC[, t, ]

    for (a in getNames(magEATdiets, dim = 2)) {
      alcohol_gap_reg <- where(magEATdiets[, t, "alcohol"][, , a] > balance_post_al[, t, "BMK"][, , a])$true$regions
      alcohol_gap_cal <- collapseNames(balance_post_al[alcohol_gap_reg, t, "BMK"][, , "2100kcal"][, , a] / dimSums(magEATdiets[alcohol_gap_reg, t, "alcohol"][, , "BMK"][, , "2100kcal"][, , a], dim = 3.4))
      alcohol_gap_cal[which(alcohol_gap_cal > 1)] <- 1

      magEATdiets[alcohol_gap_reg, t, "alcohol"][, , "BMK"][, , a] <- collapseNames(magEATdiets[alcohol_gap_reg, t, "alcohol"][, , "BMK"][, , a] * alcohol_gap_cal[alcohol_gap_reg, t, ])
    }
  }


  ### brans:
  # currently, intake of brans is set to Zero,
  # since there is no explicit information available in the EAT Lancet data set
  magEATdiets[, , "brans"] <- 0

  # Ideas for further improvement:
  # in EAT Lancet diet, cereals should be consumed in the form of wholegrain meals
  # this can be realized re-allocating a share of calories from cereals to brans according to the milling share


  ##### Mapping of EAT Lancet food commodities to MAgPIE commodities finished ####
  ################################################################################



  ##### Intermediate check if distribution from EAT Lancet food commodities (disregarding "othrcrp"and "othrcal")
  # to MAgPIE commodities is accurate

  eatLancetintakeNObalancing <- dimSums(EATdiets[, , c("othrcrp", "othrcal"), invert = TRUE],
                                        dim = 3.4)
  MAgPIEeatintakeNObalancing <- dimSums(magEATdiets[, , magSugrcrp, invert = TRUE][, , "alcohol",
                                                                                   invert = TRUE],
                                        dim = 3.4)
  checkIntakeNObalancing     <- eatLancetintakeNObalancing - MAgPIEeatintakeNObalancing

  if (max(abs(checkIntakeNObalancing)) > 1e-09) {
    warning("Distribution of EAT Lancet diet categories to MAgPIE categories below required accuracy")
  }


  ##### Calibration of mapped MAgPIE diets to total EAT Lancet Intake

  out <- magEATdiets

  if (calib) {
    eatLancetintake   <- dimSums(EATdiets, dim = 3.4)
    magEATintakeTotal <- dimSums(magEATdiets, dim = 3.4)

    intakeCalibFactor <- eatLancetintake / magEATintakeTotal
    magEATdietsCalib  <- intakeCalibFactor * magEATdiets

    out <- magEATdietsCalib
  }


  #### Set values for all countries to ZERO where no FAO statistics exist, e.i. where fsupplyHist is ZERO

  if (FAOcountr) {
    nonFAOSTAT <- where(dimSums(fsupplyHist[, "y2010", "kcal"], dim = 3.1) == 0)$true$regions
    out[nonFAOSTAT, , ] <- 0
  }


  #### Select nutrition attributes for which data should be returned

  out <- collapseNames(out[, , attributes])


  #### Define weights and units
  weight <- collapseNames(calcOutput("Population", scenario = "SSP2", aggregate = FALSE)[, years, ])
  unit   <- "kcal or wm per capita per day"

  return(list(x = out,
              weight = weight,
              unit = unit,
              description = paste0("Daily per capita food intake for MAgPIE commodities ",
                                   "consistent with EAT Lancet diet scenarios"))
  )
}
