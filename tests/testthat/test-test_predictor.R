context("tests test_predictor")

#load example data objects
#contains cm, cm_new, pheno, and regiondata
data(sysdata, package='phenopredict')

number=5
inputdata = filter_regions(expression=cm[1:10,1:10], regiondata=regiondata[1:10] ,phenodata=pheno[1:10,], phenotype="Sex", covariates=c("AGE","BMI"),type="factor", numRegions=number)

num2 = 2
predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno[1:10,], phenotype="Sex", covariates=NULL,type="factor", numRegions=num2)

predictions_test<-test_predictor(inputdata=inputdata ,phenodata=pheno[1:10,], phenotype="Sex", covariates=NULL,type="factor",predictordata=predictor)	

test_that("test_predictor() output makes sense", {
  expect_equal(length(predictions_test), 3)
  expect_equal(length(predictions_test$actual), length(predictions_test$predicted))
  expect_true(predictions_test$summarized[,"number_correct"]<=length(predictions_test$actual))
})