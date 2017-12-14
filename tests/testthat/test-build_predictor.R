context("tests build_predictor")

#load example data objects
#contains cm, cm_new, pheno, and regiondata
data(sysdata, package='phenopredict')

number=5
inputdata = filter_regions(expression=cm[1:10,1:10], regiondata=regiondata[1:10] ,phenodata=pheno[1:10,], phenotype="Sex", covariates=c("AGE","BMI"),type="factor", numRegions=number)

num2 = 2
predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno[1:10,], phenotype="Sex", covariates=NULL,type="factor", numRegions=num2)

test_that("if filter_regions() output are correct class", {
  expect_is(predictor$regiondata, 'GRanges')
  expect_is(predictor$trainingProbes, 'numeric')
  expect_is(predictor, 'list')

})

test_that("output has expected dimensions", {
  expect_equal(length(predictor), 3)
  expect_equal(length(predictor$trainingProbes), length(predictor$regiondata))
})