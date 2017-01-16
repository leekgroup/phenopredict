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
#' @param predictordata object output from build_predictor \code{predictordata}
#'
#' @return  predicted A vector of predicted phenoyptes
#'
#' @keywords phenotype, prediction, data set
#'
#' @export
#'


predict_pheno <- function(expression=NULL, phenodata=NULL, phenotype=NULL , predictordata = NULL){	
	require(minfi)
	#define possible predictions
	possibles = levels(droplevels(phenodata[,phenotype]))

	## to chose max value, but assign NA if max is 0
	which.highest <- function(x){
	    if(max(x)!=0){
	        return(which.max(x))
	    }else{
	        return(length(possibles)+1)
	    }
	}

	possNA = c(possibles,"Unassigned")
	minfi:::projectCellType(Y=expression, coefCellType=predictordata$coefEsts) -> predictions
	maxs <- apply(predictions,1,max)
	esttype = apply(predictions,1,which.highest)
	predicted <- possNA[esttype]

	return(predicted)   
}