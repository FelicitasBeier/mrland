#' @title readZabel2014
#' @description Download Zabel 2014 data
#' @author Patrick Rein
#' @importFrom utils download.file bibentry person
#' @importFrom withr local_options

downloadZabel2014 <- function() {

  local_options(timeout = 10000) # File is 7GB -> 1s/MB + buffer
  download.file("https://zenodo.org/records/3748350/files/cropsuitability.zip", "cropsuitability.zip", )
  unzip("cropsuitability.zip")
  unlink("cropsuitability.zip")

  return(list(
    title = "Global Agricultural Land Resources - A High Resolution Suitability Evaluation and Its Perspectives until 2100 under Climate Change Conditions", #nolint: line_length_linter
    url = "https://zenodo.org/records/3748350",
    doi = "10.5281/zenodo.3748350",
    version = "v2.0",
    author  = person("Florian", "Zabel"),
    license = "CC BY 4.0 Int",
    reference =
      bibentry(
        "Article",
        title = "Global Agricultural Land Resources - A High Resolution Suitability Evaluation and Its Perspectives until 2100 under Climate Change Conditions", #nolint: line_length_linter
        author = c(
          person("Florian", "Zabel"),
          person("Brigitta", "Putzenlechner"),
          person("Wolfram", "Mauser")
        ),
        year = "2014",
        journal = "PLOS ONE",
        url = "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0107522",
        doi = "10.1371/journal.pone.0107522"
      )
  ))

}
