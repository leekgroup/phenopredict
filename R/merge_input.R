
#' Uses output from select_regions() to merge output before building the 
#' predictor   
#' 
#' This function takes the output from multiple executions of select_regions() 
#' and merges the output before building and running the predictor. Objects
#' output are the merged output from select_regions but contain the same 
#' objects. 
#'
#' @param inputdata_list list out output objects from select_regions() 
#' \code{inputdata_list}
#'
#' @return inputdata merged inputdata across multiple runs of select_regions()
#'
#' @keywords phenotype, merge, prediction, expression
#'
#' @export

merge_input <- function(inputdata_list=NULL){	
	require(tidyverse)
	require(GenomicRanges)
		map(inputdata_list, function(x){return(x$regioninfo)}) %>% bind_rows -> regioninfo 
		map(inputdata_list, function(x){return(x$covmat)}) %>% bind_rows -> covmatrix
		map(inputdata_list, function(x){return(x$regions)}) %>% GRangesList %>% unlist -> regions
		
		# this also works:
		# do.call(c, map(inputdata_list, function(x){return(x$regions)}) ) -> regions2 
		
		res <- list(regioninfo = regioninfo, covmat=covmatrix, regions = regions)
		return(res)
}

