#' @title readREMIND
#' @description Reads in a reporting mif file from REMIND
#'
#' @param subtype A string composed of three items: unit, revision and indicator.
#' Unit can be either "intensive" or "extensive", revision is the input data
#' revision, and indicator is the name of the REMIND indicator
#' @return MAgPIE object with regional aggregation of REMIND H12
#' @author David Klein
#' @seealso
#' \code{\link[madrat]{readSource}}
#' @examples
#'
#' \dontrun{
#' readSource("REMIND",aggregate=FALSE)
#' }
#' @importFrom magclass read.report

readREMIND <- function(subtype) {

  .readAndRename <- function(fileList, pattern, replacement, indicator) {
    x <- NULL
    for (f in fileList) {
      # select REMIND only since newer coupled REMIND reportings also contain MAgPIE

      report <- read.report(f,
                            as.list = FALSE,
                            showSeparatorWarning = FALSE)[, , "REMIND"]

      if (indicator %in% getNames(report, dim = 3)) {
        tmp <- report[, , indicator]
      } else {
        tmp <- report[, , 1]
        tmp[] <- NA
        getNames(tmp, dim = 3) <- indicator  # assign the missing indicator name
      }

      x <- mbind(x, tmp)

    }

    # if data contains US$2005 convert them
    if (any(grepl("US\\$2005", getItems(x, dim = 3.3)))) {
      #convert from USD05MER to USD17MER based on USA values for all countries as the CO2 price is global.
      x <- x * round(GDPuc::toolConvertSingle(1, "USA", unit_in = "constant 2005 US$MER",
                                              unit_out = "constant 2017 US$MER"), 2)
    }

    # remove model and variable name
    x <- collapseNames(x)
    # shorten names of the REMIND scenarios
    getNames(x) <- gsub(pattern, replacement, getNames(x))
    return(x)
  }

  if (grepl("_", subtype)) {

    subtype  <- strsplit(subtype, split = "_")
    revision <- numeric_version(unlist(subtype)[2])
    indicator <- unlist(subtype)[3]

    # for data that was added with revisions < 4.118 always look for US$2005
    indicator <- gsub("US\\$2017", "US$2005", indicator)


    # reading in REMIND reporting files:

    # Please refer to the 2019-R2M41/readme.txt for the source of the data
    fileList <- c("2019-R2M41/REMIND_generic_r8473-trunk-C_Budg600-rem-5.mif",
                  "2019-R2M41/REMIND_generic_r8473-trunk-C_Budg950-rem-5.mif",
                  "2019-R2M41/REMIND_generic_r8473-trunk-C_Budg1300-rem-5.mif",
                  "2019-R2M41/REMIND_generic_r8473-trunk-C_NDC-rem-5.mif",
                  "2019-R2M41/REMIND_generic_r8473-trunk-C_NPi-rem-5.mif")

    out <- .readAndRename(fileList = fileList,
                          pattern = "r8473-trunk-C_",
                          replacement = "R2M41-SSP2-",
                          indicator = indicator)

    if (revision > 4.58) {
      # Please refer to the 2021-R21M42/readme.txt for the source of the data
      fileList <- c("2021-R21M42/REMIND_generic_C_SDP-NPi-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SDP-PkBudg900-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SDP-PkBudg1000-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SDP-PkBudg1100-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP1-NPi-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP1-PkBudg900-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP1-PkBudg1100-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP1-PkBudg1300-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP2-NPi-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP2-PkBudg900-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP2-PkBudg1100-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP2-PkBudg1300-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP5-NPi-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP5-PkBudg900-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP5-PkBudg1100-rem-5.mif",
                    "2021-R21M42/REMIND_generic_C_SSP5-PkBudg1300-rem-5.mif")
      if (revision > 4.103) {
        fileList <- c(fileList,
                      "2021-R21M42/REMIND_generic_C_SDP-NDC-rem-5.mif",
                      "2021-R21M42/REMIND_generic_C_SSP1-NDC-rem-5.mif",
                      "2021-R21M42/REMIND_generic_C_SSP2-NDC-rem-5.mif",
                      "2021-R21M42/REMIND_generic_C_SSP5-NDC-rem-5.mif")
      }

      out <- mbind(out, .readAndRename(fileList = fileList,
                                       pattern = "C_(SDP|SSP)",
                                       replacement = "R21M42-\\1",
                                       indicator = indicator))
    }

    if (revision >= 4.96) {
      # Please refer to the 2023-R32M46/readme.txt for the source of the data
      fileList <- c("2023-R32M46/REMIND_generic_C_SDP_MC-NDC-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SDP_MC-NPi-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SDP_MC-PkBudg650-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP1-NDC-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP1-NPi-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP1-PkBudg1050-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP1-PkBudg650-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP2EU-NDC-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP2EU-NPi-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP2EU-PkBudg1050-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP2EU-PkBudg650-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP5-NDC-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP5-NPi-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP5-PkBudg1050-rem-5.mif",
                    "2023-R32M46/REMIND_generic_C_SSP5-PkBudg650-rem-5.mif")

      out <- mbind(out, .readAndRename(fileList = fileList,
                                       pattern = "C_(SDP_MC|SSP)",
                                       replacement = "R32M46-\\1",
                                       indicator = indicator))
    }

    # for data that was added with revisions >= 4.118 look for US$2017
    indicator <- gsub("US\\$2005", "US$2017", indicator)

    if (revision >= 4.118) {
      # Please refer to the 2025-R34M410/readme.txt for the source of the data
      fileList <- c("2025-R34M410/REMIND_generic_C_SSP1-NPi2025-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP1-PkBudg1000-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP1-PkBudg650-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2_lowEn-NPi2025-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2_lowEn-PkBudg1000-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2_lowEn-PkBudg650-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2-NPi2025-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2-PkBudg1000-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP2-PkBudg650-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP3-NPi2025-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP3-PkBudg1000-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP3-rollBack-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP5-NPi2025-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP5-PkBudg1000-rawluc-rem-5.mif",
                    "2025-R34M410/REMIND_generic_C_SSP5-PkBudg650-rawluc-rem-5.mif")

      out <- mbind(out, .readAndRename(fileList = fileList,
                                       pattern = "C_(SSP)",
                                       replacement = "R34M410-\\1",
                                       indicator = indicator))
    }

    if (revision >= 4.122) {
      # Please refer to the 2025-R34withBC/readme.txt for the source of the data
      fileList <- c("2025-R34withBC/REMIND_generic_PB650-BCdef-CTS01-BM70-noFuel-CHP17.mif",
                    "2025-R34withBC/REMIND_generic_PB650-BCdef-CTS01-noFuel-CHP17.mif",
                    "2025-R34withBC/REMIND_generic_PB650-BCdef-noFuel-CHP17.mif",
                    "2025-R34withBC/REMIND_generic_PB650-BCpess-CTS01-BM70-noFuel-CHP17.mif")

      out <- mbind(out,
                   .readAndRename(fileList = fileList,
                                  pattern = "PB650-(BCdef-CTS01-BM70|BCdef-CTS01|BCdef|BCpess-CTS01-BM70)-noFuel-CHP17",
                                  replacement = "R34BC-SSP2-PkBudg650-\\1",
                                  indicator = indicator))
    }
  }

  # shorten names of the REMIND scenarios
  getNames(out) <- gsub("(-rawluc|)-rem-5", "", getNames(out))

  return(out)

}
