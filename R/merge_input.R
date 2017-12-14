
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
#'
#' @examples
#'
#' library('GenomicRanges')
#' library('dplyr')
#' ## make up some expression data for 9 rows and 30 people
#' data(sysdata, package='phenopredict')
#' ## includes R object 'cm'   
#'
#' ## Make up some some region data
#' regions <- GRanges(seqnames = 'chr2', IRanges(
#'       start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
#'       end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))
#'
#' ## make up some expression data for 9 rows and 30 people
#' exp= cm[1:length(regions),1:30]
#'
#' ## Generate first object to be merged
#' inputdata1 <- list()
#' inputdata1$covmat = exp
#' inputdata1$regiondata = regions
#'
#' ## Generate scond object to be merged
#' regions2 = GRanges(seqnames = 'chr9', IRanges(
#'       start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
#'       end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))
#' exp2 = cm[1:length(regions2),1:30]
#' inputdata2 <- list()
#' inputdata2$covmat = exp2
#' inputdata2$regiondata = regions2
#'
#' ## merge objects
#' inputdata_merged<-merge_input(inputdata_list=list(inputdata1, inputdata2))

merge_input <- function(inputdata_list=NULL){	
	requireNamespace("dplyr", quietly=TRUE)
	requireNamespace("purrr", quietly=TRUE)
	requireNamespace("plyr", quietly=TRUE)
	requireNamespace("GenomicRanges", quietly=TRUE)	
		purrr::map(inputdata_list, function(x){return(x$covmat)}) %>% plyr::ldply(., data.frame) -> covmatrix
		purrr::map(inputdata_list, function(x){return(x$regiondata)}) %>% GRangesList %>% unlist -> regiondata
		# this also works (Thanks, Leo!):
		# do.call(c, map(inputdata_list, function(x){return(x$regiondata)}) ) -> regiondata
		#regioninfo = regioninfo,
		res <- list( covmat=covmatrix, regiondata = regiondata)
		return(res)
}