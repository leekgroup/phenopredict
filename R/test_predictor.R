
#' Test accuracy of predictor on known phenotypes
#'
#' This function takes the expression data input to 
#' build_predictor() and the coefficient estimates from
#' build_predictor() for phenotype prediction. The known 
#' phenotypes are also input for comparison and 
#' asseessment of predictor accuracy.
#'
#' @param inputdata output from select_regions() \code{inputdata}
#' @param phenodata data set with phenotype information; samples in rows, 
#' variables in columns \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, binary, factor)
#' \code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param predictordata object output from build_predictor \code{predictordata}
#'
#' @return list of actual and predicted phenotype, and summarization of output
#'
#' @keywords phenotype, prediction, test
#'
#' @export
#' 

	

test_predictor <- function(inputdata=NULL ,phenodata=NULL, phenotype=NULL, covariates=NULL,type="factor",predictordata=NULL){	
	require(minfi)

	 if(is.null(inputdata)) {
	  	stop('Must specify inputdata to use. This is the output from select_regions()')
	  }
	   if(is.null(phenodata)) {
	  	stop('Must include phenotype file.')
	  }
 	if(is.null(predictordata)) {
	  	stop('Must specify predictor data to use. This is the output from build_predictor()')
	  }

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

	#extract regions
	expressiondata = inputdata$covmat[predictor$trainingProbes,]

	## predictions
	possNA = c(possibles,"Unassigned")
	minfi:::projectCellType(Y=expressiondata, coefCellType=predictordata$coefEsts) -> predictions
	maxs <- apply(predictions,1,max)
	esttype = apply(predictions,1,which.highest)
	predicted <- possNA[esttype]  
	
	#summarize data
	actual <- phenodata[,phenotype]
	number_match <- sum(predicted==actual)
	perc_correct = sum(predicted==actual)/length(actual)
	number_sites = nrow(predictordata$coefEsts)

	summarized = c(number_sites,number_match, perc_correct)
	names(summarized) <- c("sites_tested", "number_correct", "percent_correct")

	res <- list(actual = actual, predicted=predicted, summarized=summarized)

	return(res)

}