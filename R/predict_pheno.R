#' Makes predictions for selected phenotype  
#' 
#' This function outputs predictions for phenotype of 
#' interest using selected regions upon which predictor 
#' has been built
#'
#' @param inputdata_test Object containing expression data from test set and corresponding GenomicRanges object \code{expression}
#' @param phenodata data set with phenotype information; samples in rows, 
#' variables in columns \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, factor) 
#' @param covariates Which covariates to include in model \code{covariates}
#' @param predictordata object output from build_predictor \code{predictordata}
#'
#' @return A vector of predicted phenoyptes
#'
#' @keywords phenotype, prediction, data set
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
#' 	phenodata=pheno, phenotype="sex", covariates=NULL,type="factor", numRegions=2)
#' 
#' ## build phenotype predictor
#' predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno, 
#' 	phenotype="sex", covariates=NULL,type="factor", numRegions=2)
#'
#' ## determine resubstitution error
#' ## carry out prediction in training data set
#' predictions_test<-test_predictor(inputdata=inputdata ,phenodata=pheno, 
#' 	phenotype="sex", covariates=NULL,type="factor",predictordata=predictor)	
#'
#' ## generate new expressiondata set for prediction
#' exp_new= cm_new[1:length(regions),1:30]
#' ## extract test data
#' test_data<-extract_data(newexpression=exp_new, newregiondata=predictor$regiondata, 
#' 	predictordata=predictor)
#'
#' ## predict phenotype in test data
#' predictions <- predict_pheno(inputdata_test=test_data, phenodata=pheno, 
#' phenotype="sex", type="factor" ,covariates=NULL, predictordata = predictor)

predict_pheno <- function(inputdata_test=NULL, phenodata=NULL, phenotype=NULL, type=c("factor","numeric") ,covariates=NULL, predictordata = NULL){	
	requireNamespace("stats", quietly=TRUE)

	#check to makes sure same regions are included 
	if(length(inputdata_test$regiondata) != length(predictordata$regiondata)) {
		stop('The same regions used to build your predictor must be the regions used in your test set.')
	}
	if(length(GenomicRanges::findOverlaps(inputdata_test$regiondata,predictordata$regiondata))!=length(predictordata$regiondata)){
		stop('Some of the regions in your test expression set are not the same as those in expression set used to build your predictor')
	}
	if(!(identical(names(inputdata_test$regiondata),names(predictordata$regiondata)))){
		stop('Make sure your test expression set is in the same order as your prediction expression set!')
	}


	## to chose max value, but assign NA if max is 0
	which.highest <- function(x){
	    if(max(x)!=0){
	        return(which.max(x))
	    }else{
	        return(length(possibles)+1)
	    }
	}
	
	expression = inputdata_test$covmat
	regiondata = inputdata_test$regiondata

	type <- match.arg(type,c("factor", "numeric") )

	if(type=="factor"){
		requireNamespace("minfi", quietly=TRUE)
		# define possible predictions
		possibles = levels(droplevels(as.factor(phenodata[,phenotype])))
		possNA = c(possibles,"Unassigned")
		# make predictions
		minfi:::projectCellType(Y=expression, coefCellType=predictordata$coefEsts) -> predictions
		maxs <- apply(predictions,1,max)
		esttype = apply(predictions,1,which.highest)
		predicted <- possNA[esttype]
	}
	if(type=="numeric"){
		minfi:::projectCellType(Y=expression, coefCellType=predictor$coefEsts) -> predictions
		predicted = predictions
	}

	return(predicted)   
}