## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-------------
library(phenopredict)
## Track time spent on making the vignette
startTime <- Sys.time()

## ----load-packages, message = FALSE, warning = FALSE, eval = FALSE---------
#  ## load libraries
#  library('devtools')
#  install_github("leekgroup/phenopredict")
#  # document("/users/sellis/phenopredict")
#  library('phenopredict')

## ----'load-data', message = FALSE, warning = FALSE-------------------------
data(sysdata, package='phenopredict')   
#loads cm, cm_new, regiondata, pheno to run example

## ----'filter-regions', message = FALSE, warning = FALSE--------------------
# number of regions in expression data 
nrow(cm)

# number of samples included in example
ncol(cm)

inputdata<-filter_regions(expression=cm, regiondata=regiondata ,phenodata=pheno, phenotype="Sex",
    covariates=NULL,type="factor", numRegions=100)

# taking a look at output of filter_regions()
dim(inputdata$covmat)

inputdata$regiondata

head(inputdata$regioninfo)

## ----'build-predictor', message = FALSE, warning = FALSE-------------------
predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno, phenotype="Sex", covariates=NULL,type="factor", numRegions=10)

#number of probes used for prediction
length(predictor$trainingProbes)

#this contains the coefficient estimates used for prediction. 
# the number of rows corresponds to the number of sites used for prediction
# while the columns corresponds to the number of categories of your phenotype.
dim(predictor$coefEsts)


## ----'test-predictor', message = FALSE, warning = FALSE--------------------
predictions_test <-test_predictor(inputdata=inputdata ,phenodata=pheno, phenotype="Sex", 
    covariates=NULL,type="factor",predictordata=predictor)

# get summary of how prediction is doing
predictions_test$summarized

## ----'extract-data', message = FALSE, warning = FALSE----------------------
# looking at the input data for extract_data
dim(cm_new)

test_data<-extract_data(newexpression=cm_new, newregiondata=regiondata, predictordata=predictor)


## ----'predict-phenotype', message = FALSE, warning = FALSE-----------------
predictions<-predict_pheno(inputdata_test= test_data, phenodata=pheno, phenotype="Sex", covariates=NULL,type="factor", predictordata = predictor)

#looking at the output
table(predictions)

## ----reproducibility--------------------------------------------------------------------------------------------------
## Time spent creating this report:
diff(c(startTime, Sys.time()))

## Date this report was generated
message(Sys.time())

## Reproducibility info
options(width = 120)
devtools::session_info()

## ----createVignette, eval=FALSE---------------------------------------------------------------------------------------
#  ## Create the vignette
#  library('rmarkdown')
#  system.time(render('/users/sellis/phenopredict/vignettes/phenopredict-quickstart.Rmd', BiocStyle::html_document()))
#  
#  ## Extract the R code
#  library('knitr')
#  knit('/users/sellis/phenopredict/vignettes/phenopredict-quickstart.Rmd', tangle = TRUE)

