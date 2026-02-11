downloadLUH2UrbanFuture <- function() {
  .download <- function(urlArg, filename, md5) {
    if (!file.exists(filename) || tools::md5sum(filename) != md5) {
      utils::download.file(urlArg, filename, mode = "wb")
      if (tools::md5sum(filename) != md5) {
        warning("md5sum mismatch for ", filename)
      }
    }
  }

  needCertificate <- tryCatch({
    .download("https://luh.umd.edu/LUH2/LUH2_v2h/staticData_quarterdeg.nc",
              "staticData_quarterdeg.nc", "110211fdfd0f85044d203e9c71114165")
    FALSE
  }, warning = function(e) {
    TRUE
  }, error = function(e) {
    TRUE
  })

  if (needCertificate) {
    # on linux ssl verification failed with the default setup, so explicitly pass certificate here
    luhCertificateText <- paste0("-----BEGIN CERTIFICATE-----\n",
                                 "MIIHijCCBfKgAwIBAgIRAIfXgoVcGWAHRDlRHyki5MwwDQYJKoZIhvcNAQEMBQAw\n",
                                 "RDELMAkGA1UEBhMCVVMxEjAQBgNVBAoTCUludGVybmV0MjEhMB8GA1UEAxMYSW5D\n",
                                 "b21tb24gUlNBIFNlcnZlciBDQSAyMB4XDTI1MDQwODAwMDAwMFoXDTI2MDUwOTIz\n",
                                 "NTk1OVowaDELMAkGA1UEBhMCVVMxETAPBgNVBAgTCE1hcnlsYW5kMSwwKgYDVQQK\n",
                                 "EyNVbml2ZXJzaXR5IG9mIE1hcnlsYW5kLUNvbGxlZ2UgUGFyazEYMBYGA1UEAxMP\n",
                                 "Z3N3ZWIyOC51bWQuZWR1MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEA\n",
                                 "z6YE++B+cIJsG43uvtuHKO18N2Au/E/TH98I976qklXbcdaYMUyJudAgx6g7f63X\n",
                                 "0xjQqIVgfGL4GpOXYoLwLq8Hbl3Ev1OTXC6XEeT558X0uSy+o2iVGhXX2DMnUe2p\n",
                                 "0PbC1+s2XbG1osx9QjxPx/Egm2g4S+ptYdvyt/KLtuhVOJrdGOMChzHAJq/w2/yc\n",
                                 "U0vsPQ5SkA9O6iKRb4NlxqfUOyLQlnzWH7ObdTq5zIIPwU218A5k+/kgS59eZ0Np\n",
                                 "89Bq6vig/14b4yfN6jVDcsweneDWDB0I8cIqDC7Qpj2iqFD02aU09jCcIPXglFik\n",
                                 "zl59HtBpdzoLhOGLOP+Z1QIDAQABo4ID0TCCA80wHwYDVR0jBBgwFoAU70wAkqb7\n",
                                 "di5eleLJX4cbGdVN4tkwHQYDVR0OBBYEFNN34nlttcb19gQUr6FXs9UaK95tMA4G\n",
                                 "A1UdDwEB/wQEAwIFoDAMBgNVHRMBAf8EAjAAMB0GA1UdJQQWMBQGCCsGAQUFBwMB\n",
                                 "BggrBgEFBQcDAjBJBgNVHSAEQjBAMDQGCysGAQQBsjEBAgJnMCUwIwYIKwYBBQUH\n",
                                 "AgEWF2h0dHBzOi8vc2VjdGlnby5jb20vQ1BTMAgGBmeBDAECAjBABgNVHR8EOTA3\n",
                                 "MDWgM6Axhi9odHRwOi8vY3JsLnNlY3RpZ28uY29tL0luQ29tbW9uUlNBU2VydmVy\n",
                                 "Q0EyLmNybDBwBggrBgEFBQcBAQRkMGIwOwYIKwYBBQUHMAKGL2h0dHA6Ly9jcnQu\n",
                                 "c2VjdGlnby5jb20vSW5Db21tb25SU0FTZXJ2ZXJDQTIuY3J0MCMGCCsGAQUFBzAB\n",
                                 "hhdodHRwOi8vb2NzcC5zZWN0aWdvLmNvbTCCAX4GCisGAQQB1nkCBAIEggFuBIIB\n",
                                 "agFoAHUAlpdkv1VYl633Q4doNwhCd+nwOtX2pPM2bkakPw/KqcYAAAGWFXe4HgAA\n",
                                 "BAMARjBEAiBQEN1pmSAWTPQBUe00L5wjCVbQDgaZ0eTXXtftXbTNZAIgTO6QHYEN\n",
                                 "+zlijBx+jk7IeGbuMHgPBJN63+UbwFf0hmsAdwAZhtTHKKpv/roDb3gqTQGRqs4t\n",
                                 "cjEPrs5dcEEtJUzH1AAAAZYVd7fVAAAEAwBIMEYCIQC3CBzsc6UY9FLoATzQBytT\n",
                                 "0ednb4KXo7U2reRmpCPjAgIhANplSIfHeUaS3doSa8MbJBjmCJxXn47AuJ7DP8ox\n",
                                 "ZFShAHYADleUvPOuqT4zGyyZB7P3kN+bwj1xMiXdIaklrGHFTiEAAAGWFXe3tgAA\n",
                                 "BAMARzBFAiEAw7xLIO5r7NR/5gokg+wK673XXjkg2J3bplkbFb2Rls4CICOF9KcU\n",
                                 "1oxbNtloBRIj+zmpNld1esSkNhrYRHcykzSjMIHMBgNVHREEgcQwgcGCD2dzd2Vi\n",
                                 "MjgudW1kLmVkdYIYY2FyYm9ubW9uaXRvcmluZy51bWQuZWR1ggtnZWwudW1kLmVk\n",
                                 "dYISZ3N3ZWIyOHZoMS51bWQuZWR1ghJnc3dlYjI4dmgyLnVtZC5lZHWCEmdzd2Vi\n",
                                 "Mjh2aDMudW1kLmVkdYILbHVoLnVtZC5lZHWCHHd3dy5jYXJib25tb25pdG9yaW5n\n",
                                 "LnVtZC5lZHWCD3d3dy5nZWwudW1kLmVkdYIPd3d3Lmx1aC51bWQuZWR1MA0GCSqG\n",
                                 "SIb3DQEBDAUAA4IBgQAwWSeY88IDXTxRTd9geZn8oAGCoRK6WMzhNOgi4pBs3TeU\n",
                                 "X3rfp/C+QwVyxFqiL0YldAxSLvG8aKFQ04I9IdsoIOat2ubTSEouw8y7UQqhOdLS\n",
                                 "T6sEZEK9kTz+7CEt3ubEY8+N1ZVcEGmqH7taO1UsEOz96mgCh2q4iFpiRKzDH/4h\n",
                                 "yuGf4s8VcegAciyd0AN0hIrNzIW0KVSCtLtAv0wV2JkeSkFuFiaOHaZ9n4j393Uk\n",
                                 "qlg1fwm0Eu/HLFvXsotLyFMXIp1cJSWswub7PuRM5UXGEBHd356tkFMDPkYWA6Tl\n",
                                 "utFn+dAWe/wyhHu89DKF7DhH7RJmaGQypu90821ceEVCHtBUe/A0jRRURnE3BSBa\n",
                                 "LJY3+v5EPs3yLQflg30W7rJUhaYQ6cj9aqq1dUjY7iuRbAK7Ov6Q1ivcn6Q3Kxid\n",
                                 "dhFzm3bSgGUcUsxnhZ0i170IQJgC6SQ9PlCtkPT+Rmh1zrux7qn+JrNVmggtp9sU\n",
                                 "RJLkfVUTVO2PtXfRDds=\n",
                                 "-----END CERTIFICATE-----\n")
    luhCertificate <- withr::local_tempfile(fileext = ".pem")
    writeLines(luhCertificateText, luhCertificate)

    .download <- function(urlArg, filename, md5) {
      if (!file.exists(filename) || tools::md5sum(filename) != md5) {
        utils::download.file(urlArg,
                             filename,
                             mode = "wb",
                             method = "curl",
                             extra = paste0("--cacert ", luhCertificate, " --no-progress-meter"))
        if (tools::md5sum(filename) != md5) {
          warning("md5sum mismatch for ", filename)
        }
      }
    }
    .download("https://luh.umd.edu/LUH2/LUH2_v2h/staticData_quarterdeg.nc",
              "staticData_quarterdeg.nc", "110211fdfd0f85044d203e9c71114165")
  }

  getUrl <- function(model, urlPart) {
    return(paste0("https://luh.umd.edu/LUH2/LUH2_v2f/", model,
                  "/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-", urlPart))
  }

  .download(getUrl("IMAGE_SSP1_RCP19", "IMAGE-ssp119-2-1-f_gn_2015-2100.nc"),
            "states_ssp1_19.nc", "9aacafe0b06df7ead4a65dd494f1b23d")
  .download(getUrl("MESSAGE", "MESSAGE-ssp245-2-1-f_gn_2015-2100.nc"),
            "states_ssp2_45.nc", "8d2cd5ca241678c56d8b92464f79d97e")
  .download(getUrl("AIM", "AIM-ssp370-2-1-f_gn_2015-2100.nc"),
            "states_ssp3_70.nc", "18a7d1f1d7f687606fe6f7cba476d925")
  .download(getUrl("GCAM60", "GCAM-ssp460-2-1-f_gn_2015-2100.nc"),
            "states_ssp4_60.nc", "aad9baa5d03c6efc155f027387b118ca")
  .download(getUrl("MAGPIE", "MAGPIE-ssp585-2-1-f_gn_2015-2100.nc"),
            "states_ssp5_85.nc", "5f8b277e03c02fca569158a6234f2c3e")
}
