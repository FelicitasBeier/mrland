x <- new.magpie(c("A.1", "B.2"), "y2020",
                paste(rep(c("s1", "s2"), each = 4), c("crop", "past", "forest", "other"), sep = "."),
                fill = 0)
x[, , "crop"] <- 1
x[, , "past"] <- 2
x[, , "forest"] <- 10
x[, , "other"] <- 20
landUse9 <- new.magpie(c("A.1", "B.2"), "y2020",
                       c("crop", "past", "range", "forestry", "primforest", "secdforest",
                         "urban", "primother", "secdother"), fill = 0)
landUse9["A.1", , "primforest"] <- 3
landUse9["A.1", , "secdforest"] <- 1
landUse9["A.1", , "primother"] <- 1
landUse9["A.1", , "secdother"] <- 4

test_that("item order matches MAgPIE land classes", {
  seven <- toolSplitLandClasses(x, landUse9, "seven")
  expect_identical(getItems(seven, dim = 3.1), c("s1", "s2"))
  expect_identical(getItems(seven, dim = 3.2),
                   c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other"))
  nine <- toolSplitLandClasses(x, landUse9, "nine")
  expect_identical(getItems(nine, dim = 3.2),
                   c("crop", "past", "range", "forestry", "primforest", "secdforest",
                     "urban", "primother", "secdother"))
})

test_that("forest and other land are split by LanduseInitialisation shares", {
  nine <- toolSplitLandClasses(x, landUse9, "nine")
  expect_equal(as.vector(nine["A.1", , "s1.primforest"]), 7.5, tolerance = 1e-8)
  expect_equal(as.vector(nine["A.1", , "s1.secdforest"]), 2.5, tolerance = 1e-8)
  expect_equal(as.vector(nine["A.1", , "s2.primother"]), 4, tolerance = 1e-8)
  expect_equal(as.vector(nine["A.1", , "s2.secdother"]), 16, tolerance = 1e-8)
  expect_equal(as.vector(nine[, , "range"]), as.vector(x[, , "past"]))
  expect_true(all(nine[, , c("past", "forestry", "urban")] == 0))
})

test_that("all land goes to secondary where LanduseInitialisation has none", {
  nine <- toolSplitLandClasses(x, landUse9, "nine")
  expect_true(all(nine["B.2", , c("primforest", "primother")] == 0))
  expect_equal(as.vector(nine["B.2", , "secdforest"]), as.vector(x["B.2", , "forest"]))
  expect_equal(as.vector(nine["B.2", , "secdother"]), as.vector(x["B.2", , "other"]))
})

test_that("zero classes stay 0 where input is NA", {
  xNA <- x
  xNA["B.2", , ] <- NA
  for (n in c("seven", "nine")) {
    out <- toolSplitLandClasses(xNA, landUse9, n)
    expect_true(all(out[, , c("forestry", "urban")] == 0))
    expect_true(all(is.na(out["B.2", , c("crop", "secdforest")])))
  }
  expect_true(all(toolSplitLandClasses(xNA, landUse9, "nine")[, , "past"] == 0))
})

test_that("invalid nclasses fails", {
  expect_error(toolSplitLandClasses(x, landUse9, "five"), "does not exist")
})
