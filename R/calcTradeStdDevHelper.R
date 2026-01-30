#' calcTradeStdDevHelper
#'
#' Calculate rolling standard-deviation statistics for trade ratios observed over a window of years.
#'
#' Given a tidy data.frame/tibble with columns for exporter (ex), importer (im),
#' item (ItemCodeItem), year and a numeric value column (default `.value`),
#' this function computes a rolling standard deviation over `years` for each group
#' and returns the mean, max and min of that rolling standard deviation per group.
#'
#' @param dataIn data.frame or tibble; tidy table with columns specified in `group_vars`, `year_col` and `value_col`.
#' @param sdYears integer; window length (number of years) for the rolling standard deviation (must be >= 1).
#' @param groupVars character; grouping columns. Default: c("ex", "im", "ItemCodeItem").
#' @param yearCol character; name of the year column. Default: "Year".
#' @param valueCol character; name of the numeric value column. Default: ".value".
#' @return tibble with one row per group and columns: meansd<years>, maxsd<years>, minsd<years>.
#' @examples
#' \dontrun{
#' ratiodf <- as.data.frame(collapseNames(ratio), rev = 2)
#' toolTradeStdDev(ratiodf, years = 5)
#' }
#' @author David M Chen
#' @export
#' @importFrom dplyr arrange group_by mutate summarise across all_of .data
#' @importFrom zoo rollapply
calcTradeStdDevHelper <- function(dataIn,
                                  sdYears,
                                  groupVars = c("ex", "im", "ItemCodeItem"),
                                  yearCol = "Year",
                                  valueCol = ".value") {
  stopifnot(length(sdYears) == 1 && is.numeric(sdYears) && sdYears >= 1)

  meanCol <- paste0("meansd", sdYears)
  maxCol  <- paste0("maxsd", sdYears)
  minCol  <- paste0("minsd", sdYears)

  dataIn <- as.data.frame(collapseNames(dataIn), rev = 2)

  out <- dataIn %>%
    dplyr::arrange(.data[[yearCol]]) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(groupVars))) %>%
    dplyr::mutate(
      rollsd = zoo::rollapply(.data[[valueCol]],
                              width = sdYears,
                              FUN = stats::sd,
                              fill = NA,
                              align = "right")
    ) %>%
    dplyr::summarise(
      !!meanCol := mean(.data[["rollsd"]], na.rm = TRUE),
      !!maxCol  := max(.data[["rollsd"]],  na.rm = TRUE),
      !!minCol  := min(.data[["rollsd"]],  na.rm = TRUE),
      .groups = "drop"
    )

  out <- as.magpie(out,  spatial = c(1, 2), temporal = NULL)

  return(list(x = out,
              weight = NULL,
              unit = "ratio",
              isocountries = FALSE,
              description = "countries' import supply ratio. Imports/Domestic supply"))
}
