#' Makes predictions for selected phenotype  
#' 
#' This function outputs predictions for phenotype of 
#' interest using selected regions upon which predictor 
#' has been built
#'
#' @param inputdata_test Object containing expression data from test set and corresponding GRanges object \code{expression}
#' @param inputdata_test data set with phenotype information; samples in rows, 
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


predict_pheno <- function(inputdata_test=NULL, phenodata=NULL, phenotype=NULL, type=c("factor","numeric") ,covariates=NULL, predictordata = NULL){	

	#check to makes sure same regions are included 
	if(length(inputdata_test$regiondata) != length(predictordata$regiondata)) {
		stop('The same regions used to build your predictor must be the regions used in your test set.')
	}
	if(length(findOverlaps(inputdata_test$regiondata,predictordata$regiondata))!=length(predictordata$regiondata)){
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
		require(minfi)
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
		require(caret)
		phenouse <- as.numeric(phenodata[,phenotype])
		datar = as.data.frame(t(expression))
		colnames(datar) = names(regiondata)
		datar$phenouse = phenouse
		
		if(!is.null(covariates)){
			covars = phenodata[,covariates, drop=F]
			datar = cbind(datar, covars)
		}
		# predict
		plsClasses <- predict(predictordata$coefEsts, newdata = datar)
		predicted = plsClasses
	}

	return(predicted)   
}