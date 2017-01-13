context("tests select_regions")

#load example data objects
#contains cm, cm_new, pheno, and regiondata
data(sysdata, package='phenopredict')

inputdata = select_regions(expression=cm[1:10,1:10], regiondata=regiondata[1:10] ,phenodata=pheno[1:10,], phenotype="Sex", covariates=c("AGE","BMI"),type="factor", numRegions=5)

test_that("if select_regions() output are correct class", {
  expect_is(inputdata$regiondata, 'GRanges')
  expect_is(inputdata$covmat, 'matrix')
  expect_is(inputdata$regioninfo, 'data.frame')
})

test_that("output has expected dimensions", {
  expect_equal(length(inputdata$regiondata), nrow(inputdata$covmat))
  expect_equal(nrow(inputdata$covmat), nrow(inputdata$regioninfo))
})