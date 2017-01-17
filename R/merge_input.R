
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
#' @return merged inputdata across multiple runs of select_regions()
#'
#' @keywords phenotype, merge, prediction, expression
#'
#' @export

merge_input <- function(inputdata_list=NULL){	
	require(tidyverse)
	require(GenomicRanges)	
		map(inputdata_list, function(x){return(x$covmat)}) %>% ldply(., data.frame) -> covmatrix
		map(inputdata_list, function(x){return(x$regioninfo)}) %>% ldply(., data.frame) -> regioninfo
		# map(inputdata_list, function(x){return(x$regioninfo)}) %>% bind_rows -> regioninfo 
		# map(inputdata_list, function(x){return(x$covmat)}) %>% bind_rows(.id=NULL) -> covmatrix
		map(inputdata_list, function(x){return(x$regiondata)}) %>% GRangesList %>% unlist -> regiondata
		
		# this also works:
		# do.call(c, map(inputdata_list, function(x){return(x$regiondata)}) ) -> regiondata2 
		
		res <- list(regioninfo = regioninfo, covmat=covmatrix, regiondata = regiondata)
		return(res)
}