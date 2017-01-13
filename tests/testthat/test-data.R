context("tests example data")

#load example data objects
#contains cm, cm_new, pheno, and regiondata
data(sysdata, package='phenopredict')

test_that("if example data are correct class", {
  expect_is(regiondata, 'GRanges')
  expect_is(cm, 'matrix')
  expect_is(cm_new, 'matrix')
  expect_is(pheno, 'data.frame')
})


test_that("data have expected dimensions", {
  expect_equal(length(regiondata), nrow(cm))
  expect_equal(nrow(cm), nrow(cm_new))
  expect_equal(ncol(cm), nrow(pheno))
})