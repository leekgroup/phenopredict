#' Makes predictions for selected phenotype  
#' 
#' This function outputs predictions for phenotype of 
#' interest using selected regions upon which predictor 
#' has been built
#'
#' @param expression Expression data from test set \code{expression}
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


predict_pheno <- function(expression=NULL, phenodata=NULL, phenotype=NULL, type=c("factor","numeric") ,covariates=NULL, predictordata = NULL){	
	require(minfi)


	## to chose max value, but assign NA if max is 0
	which.highest <- function(x){
	    if(max(x)!=0){
	        return(which.max(x))
	    }else{
	        return(length(possibles)+1)
	    }
	}
	
	type <- match.arg(type,c("factor", "numeric") )


	if(type=="factor"){
		# define possible predictions
		possibles = levels(droplevels(phenodata[,phenotype]))
		possNA = c(possibles,"Unassigned")
		# make predictions
		minfi:::projectCellType(Y=expression, coefCellType=predictordata$coefEsts) -> predictions
		maxs <- apply(predictions,1,max)
		esttype = apply(predictions,1,which.highest)
		predicted <- possNA[esttype]
	}
	if(type=="numeric"){
		phenouse <- as.numeric(phenodata[,phenotype])
		datar = as.data.frame(t(expression))
		datar$phenouse = phenouse
		covars = phenodata[,covariates, drop=F]
		datar = cbind(datar, covars)

		# predict
		plsClasses <- predict(predictordata$coefEsts, newdata = datar)
		predicted = plsClasses
	}

	return(predicted)   
}