#' Extracts data for predictions using regions from select_region()  
#' 
#' This function extracts regions of interest from data set on which you want 
#' to make predictions.
#'
#' @param expression expression data set from which to extract coverage data 
#' \code{expression}
#' @param inputdata object output from select_regions() \code{inputdata}
#' @param predictordata object output from build_predictor \code{predictordata}
#'
#' @return  comat_test An n x m data.frame of the selected regions from the 
#' data set specified by `expression`
#'
#' @keywords phenotype, prediction, data set
#'
#' @export
#' 


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

	  #define function to extract regions when more than one input file supplied
	  getexp <- function(x,y){
	  				out<-x[y$regioninfo$index,]
	  				return(out)
	  			}
	#if we're extracting over multiple input objects (aka if you had to use merge_input(), look here
	  if(is.list(expression)){
	  		if(!is.list(inputdata)) {
	  			stop('If >1 expression file is input, inputdata must also be a list.')
	  		}
	  		#extract regions from each expression data set before merging [using inputdata$regioninfo$index]
	  		Map(getexp,expression,inputdata) -> tmp
	  		#combine across dataframes
	  		ldply(tmp, data.frame) -> covmat_test
	#but if you only have one inputdata, it's pretty straightforward, and look here
	 }else{
		yGene=expression

		#extract appropriate regions
		covmat_test<-yGene[unique(inputdata$regioninfo$index),]
		
	}
	#either way, let's just pull out the regions we actually want to use to build the predictor
	covmat_test = covmat_test[predictor$trainingProbes,]
	return(covmat_test)
}




