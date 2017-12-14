context("tests predict_pheno")

#load example data objects
#contains cm, cm_new, pheno, and regiondata
data(sysdata, package='phenopredict')

number=5
inputdata = filter_regions(expression=cm[1:10,1:10], regiondata=regiondata[1:10] ,phenodata=pheno[1:10,], phenotype="Sex", covariates=c("AGE","BMI"),type="factor", numRegions=number)

num2 = 2
predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno[1:10,], phenotype="Sex", covariates=NULL,type="factor", numRegions=num2)

predictions_test<-test_predictor(inputdata=inputdata ,phenodata=pheno[1:10,], phenotype="Sex", covariates=NULL,type="factor",predictordata=predictor)	

exp_new= cm_new[1:10,1:10]
## extract test data
test_data<-extract_data(newexpression=exp_new, newregiondata=predictor$regiondata, 
	predictordata=predictor)

## predict phenotype in test data
predictions <- predict_pheno(inputdata_test=test_data, phenodata=pheno[1:10,], phenotype="Sex", type="factor" ,covariates=NULL, predictordata = predictor)

test_that("predict_pheno() output make sense", {
  expect_equal(length(predictions), ncol(test_data$covmat))
})