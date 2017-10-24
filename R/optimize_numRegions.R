#' Optimize number of regions used for prediction 
#'
#' This function takes a list of possible options for 
#' the numRegions argument of build_predictor. 
#' Using this set of possible numRegions and expression data 
#' (the training data / output of filter_regions()), 
#' this function builds a predictor for each possible 
#' numRegions. Prediction accuracy is then calculated  
#' across varying numbers of regions. The numRegions argument
#' that optimizes accuracy in the training data can then be 
#' used in build_predictor.
#'
#' @param inputdata output from filter_regions() \code{inputdata}
#' @param phenodata data set with phenotype information; samples in rows, 
#' variables in columns \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, binary, factor)
#' \code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param numRegions_set set of numRegions to test \code{numRegions_set}
#'
#' @return Prediction accuracies across each numRegions argument tested
#'
#' @keywords phenotype, prediction, optimization
#'
#' @export
#' 
#' @examples
#'
#'

optimize_numRegions <- function(inputdata=NULL ,phenodata=NULL, phenotype=NULL, covariates=NULL,type=NULL,numRegions_set=c(10,20,30,60,80,100,150,200)){

accuracies <- c()

for(regions in numRegions_set){
	predictor<-build_predictor(inputdata=inputdata ,phenodata=phenodata, phenotype=phenotype, covariates=NULL,type=type, numRegions=regions)

	predictions_test <-test_predictor(inputdata=inputdata ,phenodata=pheno, phenotype=phenotype,covariates=NULL,type=type,predictordata=predictor )
	if(type=="factor"){
		accuracies <- c(accuracies,as.numeric(predictions_test$summarized[,"percent_correct"]))
	}
	if(type=="numeric"){
		accuracies <- c(accuracies,as.numeric(predictions_test$summarized[,"mean_diff"]))
	}
}
names(accuracies) <- numRegions_set

if(type=="factor"){
	which_to_use<-min(which(accuracies==max(accuracies)))
	numRegions<-numRegions_set[which_to_use]
}
if(type=="numeric"){
	which_to_use<-min(which(accuracies==min(accuracies)))
	numRegions<-numRegions_set[which_to_use]
}

out = list(accuracies=accuracies, numRegions=numRegions)
return(out)
}