#' Extracts data for predictions using regions from select_region()  
#' 
#' This function extracts regions of interest from data set on which you want to make predictions.
#'
#' @param expression expression data set from which to extract coverage data \code{expression}
#' @param inputdata object output from select_regions() \code{inputdata}
#' @param predictordata object output from build_predictor \code{predictordata}
#'
#' @return  comat_test An n x m data.frame of the selected regions from the data set specified by `expression`
#'
#' @keywords phenotype, prediction, data set
#'
#' @export
#' 
#' @examples
#' testdata <-extract_data(expression=cm_new, inputdata=inputdata, predictordata=predictor)

extract_data <- function(expression=NULL, inputdata=NULL, predictordata=NULL){	
	##########
	### Get the unique regions
	### from the GTEX selection proces
	##########
	## first, some checks
		##ADD thes
	 if(is.null(expression)) {
	  	stop('Expression Data must be supplied.')
	  }
	   if(is.null(inputdata)) {
	  	stop('Must specify inputdata to use. This is the output from select_regions()')
	  }
 	if(is.null(predictordata)) {
	  	stop('Must specify predictor data to use. This is the output from build_predictor()')
	  }

	yGene=expression

	#extract appropriate regions
	covmat_test<-yGene[unique(inputdata$regioninfo$index),]
	covmat_test = covmat_test[predictor$trainingProbes,]
	
	return(covmat_test)
}