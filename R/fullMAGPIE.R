#' fullMAgPIE
#'
#' Function that produces the regional data set for running the MAgPIE model.
#'
#' @param rev data revision which should be used as input (numeric_version).
#' @param dev For developing purposes, apply changes as per dev flag
#' @author Jan Philipp Dietrich, Benjamin Leon Bodirsky, Florian Humpenoeder, Edna J. Molina Bacca
#' @seealso
#' \code{\link{readSource}}, \code{\link{getCalculations}}, \code{\link{calcOutput}}
#' @examples
#' \dontrun{
#' retrieveData("MAGPIE", rev = numeric_version("12"),
#'              mainfolder = "pathtowhereallfilesarestored")
#' }
#'
fullMAGPIE <- function(rev = numeric_version("0.1"), dev = "") {

  if (rev < numeric_version("4.66")) {
    stop("This version of mrland does not support revision below 4.66 anymore. ",
         "Please use a older snapshot/version of the library, if you need older revisions.")
  }

  magYears <- findset("time")
  magYearsPast <- findset("past")
  shortYears <- findset("t_all")

  # Check if mapping comes with an additional "superregion" layer and if so,
  # aggregate some outputs to the superregional level. Otherwise,
  # aggregate everything to regional level.
  map <- toolGetMapping(getConfig("regionmapping"), type = "regional", where = "mappingfolder")
  superregion <- ifelse("superregion" %in% colnames(map), "superregion", "region")

  cellsregions <- function(reg_revision, map) { # nolint
    # function which calculates the name vector for spatial 0.5 degree MAgPIE data sets
    # containing MAgPIE cell number and corresponding region
    regionscode    <- regionscode(map)
    spatial_header <- spatial_header(map) # nolint
    save(spatial_header, regionscode, map, reg_revision, file = "spatial_header.rda", compress = "xz")
  }
  cellsregions(rev, map)

  # data fully agrees with the data currently used in MAgPIE and new data set is implemented
  calcOutput("TauTotal",  years = 1995,        round = 2, file = "fm_tau1995.cs4", aggregate = superregion)

  # 09 drivers
  calcOutput("GDP",
             unit = "constant 2017 Int$PPP",
             naming = "scenario",
             aggregate = FALSE,
             years = magYears,
             round = 3, # please dont increase rounding, this can create errors
             file = "f09_gdp_ppp_iso.csv")
  calcOutput("GDP",
             unit = "constant 2017 US$MER",
             naming = "scenario",
             aggregate = FALSE,
             years = magYears,
             round = 3, # please dont increase rounding, this can create errors
             file = "f09_gdp_mer_iso.csv")
  calcOutput("Population",
             naming = "scenario",
             aggregate = FALSE,
             years = magYears,
             round = 6, # please dont increase rounding, this can create errors
             file = "f09_pop_iso.csv")
  calcOutput("Urban",
             asShare = FALSE,
             naming = "scenario",
             aggregate = FALSE,
             years = magYears,
             round = 4, # please dont increase rounding, this can create errors
             file = "f09_urban_iso.csv")

  calcOutput("DevelopmentState", round = 4, file = "f09_development_state.cs3")
  calcOutput("GovernanceIndicator", years = shortYears, round = 4, file = "f09_governance_indicator.cs3")
  calcOutput("Demography", education = FALSE, aggregate = FALSE, file = "f09_demography.cs3",
             round = 6) # please dont increase rounding, this can create errors
  calcOutput("PhysicalInactivity", aggregate = FALSE, years = magYears, round = 3,
             file = "f09_physical_inactivity.cs3")
  calcOutput("GDPdeflator", aggregate = FALSE, round = 3, file = "fm_gdp_defl_ppp.cs4",
             currency = "PPP")

  # 13 tc
  calcOutput("ExoTcDummy",       round = 4, file = "f13_tau_scenario.csv", aggregate = superregion)
  calcOutput("TCguess",          round = 3, file = "f13_tcguess.cs4", aggregate = superregion)
  calcOutput("TauHistorical",    round = 2, file = "f13_tau_historical.csv", aggregate = superregion)

  # 14 yields
  calcOutput("PastureYield", round = 3, file = "f14_pasture_yields_hist.csv")
  calcOutput("FAOYield", cut = 0.98, years = magYearsPast, round = 2, file = "f14_region_yields.cs3")
  calcOutput("Ir2RfYieldRatio", round = 2, file = "f14_ir2rf_ratio.cs4")

  # 15 food
  calcOutput("BodyHeight", aggregate = FALSE, years = magYearsPast, round = 2, file = "f15_bodyheight_historical.cs3")
  calcOutput("RegressionParameters", aggregate = FALSE, round = 3, file = "f15_schofield_parameters.cs3",
             regression = "Schofield")
  calcOutput("RegressionParameters", aggregate = FALSE, round = 3, file = "f15_schofield_parameters_height.cs3",
             regression = "FAO_WHO_UNU1985")
  calcOutput("RegressionParameters", aggregate = FALSE, round = 3, file = "f15_bmi_shr_regr_paras.cs3",
             regression = "bmi_shr")
  calcOutput("RegressionParameters", aggregate = FALSE, round = 3, file = "f15_demand_regression_parameters.cs3",
             regression = "demand_regression")
  calcOutput("RegressionParameters", aggregate = FALSE, round = 3, file = "f15_bodyheight_regr_paras.cs3",
             regression = "bodyheight_regression")

  calcOutput("Intake", modelinput = "age_groups_hist", standardize = FALSE, method = "FAO_WHO_UNU1985",
             aggregate = FALSE, years = magYearsPast, round = 1, file = "f15_intake_pc_observed_iso.cs3")
  calcOutput("FoodSupplyPast", per_capita = TRUE, products = NULL, product_aggr = FALSE,
             populationweight = "PopulationPast", attributes = "kcal", aggregate = FALSE,
             years = magYearsPast, round = 1, file = "f15_kcal_pc_iso.csv")
  calcOutput("Household_balanceflow",    years = magYears, round = 4, file = "f15_household_balanceflow.cs3")
  calcOutput("NutritionAttributes",      years = magYears, round = 4, file = "f15_nutrition_attributes.cs3",
             aggregate = FALSE)
  calcOutput("IniFoodPrice", datasource = "FAO", years = NULL, round = 4,
             products = "kall", file = "f15_prices_initial.csv", aggregate = FALSE, year = "y2005")
  calcOutput("BMIshr", convert = TRUE, years = magYearsPast, round = 4, file = "f15_bmi_shr_past.cs3",
             aggregate = FALSE)
  calcOutput("BMI", file = "f15_bmi.cs3", aggregate = FALSE)

  calcOutput("EATLancetDiets",  aggregate = TRUE, round = 4, file = "f15_intake_EATLancet.cs3",
             attributes = "kcal", calib = TRUE, FAOcountr = FALSE)
  calcOutput("EATLancetDiets",  aggregate = FALSE, round = 1, file = "f15_intake_EATLancet_iso.cs3",
             attributes = "kcal", calib = TRUE, FAOcountr = FALSE)
  calcOutput("EATLancetWaste",  aggregate = TRUE, round = 4, file = "f15_supply2intake_ratio_bottomup.cs3",
             out_type = "ratio_detailed")
  calcOutput("EATLancetWaste",  aggregate = FALSE, round = 4, file = "f15_supply2intake_ratio_bottomup_iso.cs3",
             out_type = "ratio_detailed")
  calcOutput("EATLancetWaste",  aggregate = TRUE, round = 4, file = "f15_calib_factor_FAOfsupply.cs4",
             out_type = "calib")
  calcOutput("EATLancetWaste",  aggregate = FALSE, round = 4, file = "f15_calib_factor_FAOfsupply_iso.cs4",
             out_type = "calib")
  calcOutput("FAOLossesWaste",  aggregate = TRUE, round = 4, file = "f15_supply2intake_ratio_FAO.cs3",
             out_type = "waste")
  calcOutput("FAOLossesWaste",  aggregate = FALSE, round = 4, file = "f15_supply2intake_ratio_FAO_iso.cs3",
             out_type = "waste")
  calcOutput("EATLancetTargets", aggregate = TRUE, round = 4, file = "f15_targets_EATLancet.cs3",
             attributes = "kcal/d")
  calcOutput("EATLancetTargets", aggregate = FALSE, round = 4, file = "f15_targets_EATLancet_iso.cs3",
             attributes = "kcal/d")
  calcOutput("EATFruitvegRatio", aggregate = TRUE, round = 4, file = "f15_fruitveg2others_kcal_ratio.cs3",
             populationweight = "PopulationPast")
  calcOutput("EATFruitvegRatio", aggregate = FALSE, round = 4, file = "f15_fruitveg2others_kcal_ratio_iso.cs3",
             populationweight = "PopulationPast")

  calcOutput("NINDiets",  aggregate = TRUE, round = 4, file = "f15_intake_NIN.cs3",
             attributes = "kcal", calib = TRUE, FAOcountr = FALSE)
  calcOutput("NINDiets",  aggregate = FALSE, round = 1, file = "f15_intake_NIN_iso.cs3",
             attributes = "kcal", calib = TRUE, FAOcountr = FALSE)

  # 16 demand
  calcOutput("Attributes", round = 4, aggregate = FALSE,        file = "fm_attributes.cs3")
  calcOutput("SeedShare", years = magYears, round = 4,           file = "f16_seed_shr.csv")
  calcOutput("LossShare", years = magYears, round = 4,           file = "f16_waste_shr.csv")
  calcOutput("DomesticBalanceflow", years = magYears, round = 4, file = "f16_domestic_balanceflow.csv")

  # 18 residues
  calcOutput("Multicropping", extend_future = TRUE, years = magYears, round = 4,
             file = "f18_multicropping.csv", aggregate = TRUE)
  calcOutput("ResCombustEff", round = 4, file = "f18_res_combust_eff.cs4", aggregate = FALSE)

  # 20 processing
  calcOutput("Processing_shares",             years = magYears, round = 4,
             file = "f20_processing_shares.cs3")
  calcOutput("Processing_conversion_factors", years = magYears, round = 4,
             file = "f20_processing_conversion_factors.cs3", aggregate = FALSE)
  calcOutput("Processing_balanceflow",        years = magYears, round = 4,
             file = "f20_processing_balanceflow.cs3")

  # 21 trade
  calcOutput("TradeSelfSuff",    years = magYears, round = 2, file = "f21_trade_self_suff.cs3",
             aggregate = superregion)
  calcOutput("TradeExportShr",   years = magYears, round = 2, file = "f21_trade_export_share.cs3",
             aggregate = superregion)
  calcOutput("TradeBalanceflow", years = magYears, round = 4, file = "f21_trade_balanceflow.cs3",
             aggregate = FALSE)
  calcOutput("TradeBalance", years = magYears, round = 2, file = "f21_trade_balance.cs3",
             aggregate = superregion)
  calcOutput("TradeMargin", round = 4, file = "f21_trade_margin.cs3",
             aggregate = superregion)
  calcOutput("TradeTariff", round = 4, file = "f21_trade_tariff.cs3",
             aggregate = superregion)
  calcOutput("TradeMargin", bilateral = TRUE, round = 4, file = "f21_trade_margin_bilat.cs5",
             aggregate = TRUE)
  calcOutput("TradeTariff", bilateral = TRUE, round = 4, file = "f21_trade_tariff_bilat.cs5",
             aggregate = TRUE)

  # 31 Past
  calcOutput("PastureYield", range_pastr = TRUE, round = 3, file = "f31_grassl_yld_hist.cs3")

  # 32 forestry
  calcOutput("AfforestCosts", years = 2001,        round = 0, file = "f32_fac_req_ha.csv")
  calcOutput("GrowingStockPlantations", aggregate = TRUE, round = 0, file = "f32_gs_target.cs4")
  calcOutput("GrowingStockPlantAbsolute", aggregate = TRUE, round = 0, file = "f32_gs_absolutetarget.cs4")
  calcOutput("GrowingStockpha", aggregate = TRUE, round = 0, file = "f32_gs_relativetarget.cs4")
  calcOutput("PlantationContribution", aggregate = TRUE, round = 3, file = "f32_plantation_contribution.cs3")
  calcOutput("PlantedForest", aggregate = TRUE, round = 3, file = "f32_plantedforest.cs4")
  calcOutput("PlantEstablishCalib", aggregate = TRUE, round = 2, file = "f32_estb_calib.cs4")
  calcOutput("TradeSelfSuff",    years = magYears, round = 2, file = "f32_trade_self_suff.cs3",
             aggregate = superregion)

  # 35 natural vegetation
  calcOutput("ForestLossShare", round = 7, file = "f35_forest_lost_share.cs3")
  calcOutput("ForestDisturbances", round = 7, file = "f35_forest_disturbance_share.cs4")
  calcOutput("GrowingStockNatVegAbsolute", aggregate = TRUE, round = 0, file = "f35_gs_absolutetarget.cs4")
  calcOutput("GrowingStockNRF", aggregate = TRUE, round = 0, file = "f35_gs_relativetarget.cs4")

  # 36 employment
  calcOutput("WeeklyHoursILO", projections = TRUE, aggregate = TRUE, years = seq(1965, 2150, 5),
             round = 2, file = "f36_weekly_hours.csv")
  calcOutput("WeeklyHoursILO", projections = TRUE, aggregate = FALSE, years = seq(1965, 2150, 5),
             round = 2, file = "f36_weekly_hours_iso.csv")
  calcOutput("HourlyLaborCosts", projection = "SSP2", aggregate = FALSE, years = seq(1965, 2015, 5),
             round = 4, file = "f36_historic_hourly_labor_costs.csv")
  calcOutput("RegressionsILO", subtype = "HourlyLaborCosts", recalculate = FALSE, aggregate = FALSE,
             round = 10, file = "f36_regression_hourly_labor_costs.csv")
  calcOutput("AgEmplILO", aggregate = FALSE, subsectors = FALSE, years = seq(1995, 2015, 5),
             round = 4, file = "f36_historic_ag_employment.csv")
  calcOutput("NonMAgPIEFactorCosts", subtype = "subsidies", aggSubsidies = TRUE, years = seq(1965, 2150, 5),
             round = 4, file = "f36_unspecified_subsidies.csv")
  calcOutput("NonMAgPIEFactorCosts", subtype = "missingVoP", years = seq(1965, 2150, 5),
             round = 4, file = "f36_nonmagpie_factor_costs.csv")

  # 38 factor costs
  calcOutput("FAOYield", cut = 0.98, years = 1995, round = 2, file = "f38_region_yield.csv")
  # Question: Is f38_region_yield used? And why not f14_region_yield?
  calcOutput("FacReq", round = 2, aggregate = "GLO", years = 2005, file = "f38_fac_req_fao.csv")
  calcOutput("FacReq", round = 2, aggregate = TRUE, file = "f38_fac_req_fao_regional.cs4")
  calcOutput("AgCapLabourShare", round = 4, aggregate = FALSE, years = c(1995, 2000, 2005, 2010),
             file = "f38_historical_share_iso.csv")
  calcOutput("RegFactorShare", datasource = "USDA", round = 4, aggregate = FALSE,
             file = "f38_regression_cap_share.csv")
  calcOutput("FactorCosts", aggregate = FALSE, years = c(1995, 2000, 2005, 2010), round = 2,
             file = "f38_hist_factor_costs_iso.csv")

  # 41 Area Equipped for Irrigation
  # f41_irrig(j) should be read out of calcAreaEquippedForIrrigation()
  calcOutput("IrrigationInvCosts", years = shortYears, round = 0, file = "f41_c_irrig.csv")

  # 42_water_demand
  calcOutput("PumpingCosts", round = 4, file = "f42_pumping_cost.cs4", aggregate = TRUE)

  # 50 n soil budget
  calcOutput("SNUpE", years = magYears, round = 4, file = "f50_snupe.cs4", rev = rev, maccbase = FALSE)
  calcOutput("SNUpE", years = magYears, round = 4, file = "f50_snupe_base.cs4", rev = rev, maccbase = TRUE)
  calcOutput("NitrogenBudgetBalanceflow", years = magYears, round = 4, file = "f50_nitrogen_balanceflow.cs4")
  calcOutput("NitrogenFixationNdfa", years = magYears, round = 4, file = "f50_ndfa.cs4")
  calcOutput("NitrogenFixationFreeliving", round = 4, file = "f50_fixation_freeliving.cs4", aggregate = FALSE)
  calcOutput("AtmosphericDepositionRates", round = 4, file = "f50_atmospheric_deposition_rates.cs4")

  calcOutput("NuePasture", years = magYears, round = 4, file = "f50_nue_pasture.cs4", maccbase = FALSE)
  calcOutput("NuePasture", years = magYears, round = 4, file = "f50_nue_base_pasture.cs4", maccbase = TRUE)
  calcOutput("NitrogenBudgetPastureBalanceflow", years = magYears, round = 4,
             file = "f50_nitrogen_balanceflow_pasture.cs4")
  calcOutput("NitrogenFixationRatePasture", years = magYears, round = 5,
             file = "f50_nitrogen_fixation_rates_pasture.cs4")

  # 51 nitrogen pollution
  calcOutput("EfNSoil",  round = 4, file = "f51_ef_n_soil.cs3", aggregate = FALSE, method = "IPCC")
  calcOutput("EfNSoil",  round = 4, file = "f51_ef_n_soil_reg.cs3")
  calcOutput("EF3confinement", round = 4, file = "f51_ef3_confinement.cs4")
  calcOutput("EF3prp", round = 4, file = "f51_ef3_prp.cs4")

  # 52 carbon
  calcOutput("AdjustGrassi2021", aggregate = TRUE, file = "f52_land_carbon_sink_adjust_grassi.cs3")

  # 53 methane
  calcOutput("EFch4Rice", years = magYears, round = 4, file = "f53_EFch4Rice.cs4")
  calcOutput("EFch4AWMS", years = magYears, round = 4, file = "f53_EFch4AWMS.cs4")

  # 55 awms
  calcOutput("ManureFuelShr",      years = magYears,  round = 4, file = "f55_manure_fuel_shr.cs4")
  calcOutput("AWMSconfShr",        years = magYears,  round = 4, file = "f55_awms_shr.cs4", rev = rev)
  calcOutput("EF3confinement", selection = "recycling", round = 4, file = "f55_awms_recycling_share.cs4")

  # 56_ghg_policy
  calcOutput("GHGPrices", datasource = "SSP_and_REM", years = shortYears, round = 2,
             file = "f56_pollutant_prices.cs3", rev = rev)
  calcOutput("GHGPrices", datasource = "S4N_project", round = 2,
             file = "f56_pollutant_prices_sim4nexus.cs3", rev = rev)

  # 57 maccs
  shortYearsFrom2010 <- shortYears[as.integer(sub("y", "", shortYears)) >= 2010]
  calcOutput("MACCsN2O", sector = "landuse", source = "ImageMacc", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_n2o.cs3")
  calcOutput("MACCsCH4", sector = "landuse", source = "ImageMacc", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_ch4.cs3")
  calcOutput("MACCsN2O", sector = "landuse", source = "PBL_MACC_2019", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_n2o_2019.cs3")
  calcOutput("MACCsCH4", sector = "landuse", source = "PBL_MACC_2019", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_ch4_2019.cs3")
  calcOutput("MACCsN2O", sector = "landuse", source = "PBL_MACC_2022", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_n2o_2022.cs3")
  calcOutput("MACCsCH4", sector = "landuse", source = "PBL_MACC_2022", years = shortYearsFrom2010,
             round = 4, file = "f57_maccs_ch4_2022.cs3")


  # 59 som
  calcOutput("SOMexogenous", years = magYears, round = 4, file = "f59_som_exogenous.cs3")
  calcOutput("SoilStockChangeFactors", round = 2, file = "f59_ch5_F_LU_2019reg.cs3")

  # 60 bioenergy
  calcOutput("1stBioDem", years = magYears, round = 3, file = "f60_1stgen_bioenergy_dem.cs3")
  calcOutput("2ndBioDem", datasource = "SSP_and_REM", years = shortYears, round = 3,
             file = "f60_bioenergy_dem.cs3", rev = rev)
  calcOutput("2ndBioDem", datasource = "S4N_project", round = 3,
             file = "f60_bioenergy_dem_sim4nexus.cs3", rev = rev)
  calcOutput("ResFor2ndBioengery", products = "kres", product_aggr = TRUE, add_off = TRUE,
             years = magYears, round = 3, file = "f60_2ndgenBE_residue_dem.cs3")

  # 62 Material
  calcOutput("DemMaterial", years = magYearsPast, round = 4, file = "f62_dem_material.cs3")
  calcOutput("BioplasticToBiomass", aggregate = FALSE, round = 4, file = "f62_bioplastic2biomass.csv")
  calcOutput("HistBioplasticProd", aggregate = FALSE, years = c(2010, 2015, 2020), round = 3,
             file = "f62_hist_dem_bioplastic.csv")

  # 70 livestock
  calcOutput("FeedBaskets",           years = magYears,  round = 4, file = "f70_feed_baskets.cs3")
  calcOutput("FeedBalanceflow",       years = magYears,  round = 4, file = "f70_feed_balanceflow.cs3")
  calcOutput("LivestockProductivity", years = magYears,  round = 4, file = "f70_livestock_productivity.cs3")
  calcOutput("SlaughterFeedShare",    years = magYears,  round = 4, file = "f70_slaughter_feed_share.cs4")
  calcOutput("PYieldSlope",                             round = 2, file = "f70_pyld_slope_reg.cs4")
  calcOutput("Production", round = 4, products = "kli", aggregate = TRUE, years = seq(1995, 2010, 5),
             file = "f70_hist_prod_livst.cs3")
  calcOutput("FactorCostsLivst", round = 4, aggregate = TRUE, years = seq(1995, 2015, 5),
             file = "f70_hist_factor_costs_livst.cs3")

  # 73 timber -- Always needed on iso country level so no need to aggregate
  calcOutput("EndUseTimber", aggregate = FALSE, round = 4, file = "f73_prod_specific_timber.csv")
  calcOutput("EndUseTimber", aggregate = TRUE, round = 2, file = "f73_regional_timber_demand.csv")
  calcOutput("ConstructionWoodDemand", aggregate = TRUE, round = 2, file = "f73_construction_wood_demand.cs3")
}
