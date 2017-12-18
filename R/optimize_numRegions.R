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
#' library('GenomicRanges')
#' library('dplyr')
#'
#' ## Make up some some region data
#' regions <- GRanges(seqnames = 'chr2', IRanges(
#'      start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
#'      end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))
#'
#' ## make up some expression data for 9 rows and 30 people
#' data(sysdata, package='phenopredict')   
#' ## includes R object 'cm'
#' exp= cm[1:length(regions),1:30]
#'
#' ## generate some phenotype information
#' sex = as.data.frame(rep(c("male","female"),each=15))
#' age = as.data.frame(sample(1:100,30))
#' pheno = dplyr::bind_cols(sex,age)
#' colnames(pheno) <- c("sex","age")
#'
#' ## filter regions to be used to build the predictor
#' inputdata <- filter_regions(expression=exp, regiondata=regions,
#' 	phenodata=pheno, phenotype="sex", covariates=NULL,type="factor", numRegions=5)
#' 
#' regnum <- optimize_numRegions(inputdata=inputdata ,phenodata=pheno, 
#' phenotype="sex", covariates=NULL,type="factor",numRegions_set=c(3,5))

optimize_numRegions <- function(inputdata=NULL ,phenodata=NULL, phenotype=NULL, 
	covariates=NULL,type=NULL,numRegions_set=c(10,20,30,60,80,100,150,200)){

accuracies <- c()


for(regions in numRegions_set){
	predictor<-build_predictor(inputdata=inputdata ,phenodata=phenodata, 
		phenotype=phenotype, covariates=NULL,type=type, numRegions=regions)

	predictions_test <-test_predictor(inputdata=inputdata,
		phenodata=phenodata, phenotype=phenotype,covariates=NULL,
		type=type,predictordata=predictor )
	if(type=="factor"){
		accuracies <- c(accuracies,as.numeric(
			predictions_test$summarized[,"percent_correct"]))
	}
	if(type=="numeric"){
		accuracies <- c(accuracies,as.numeric(
			predictions_test$summarized[,"correlation"]))
	}
}
names(accuracies) <- numRegions_set

if(type=="factor"){
	which_to_use<-min(which(accuracies==max(accuracies)))
	numRegions<-numRegions_set[which_to_use]
}
if(type=="numeric"){
	which_to_use<-min(which(accuracies==max(accuracies)))
	numRegions<-numRegions_set[which_to_use]
}

out = list(accuracies=accuracies, numRegions=numRegions)
return(out)
}