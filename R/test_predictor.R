
#' Test accuracy of predictor on known phenotypes
#'
#' This function takes the expression data input to 
#' build_predictor() and the coefficient estimates from
#' build_predictor() for phenotype prediction. The known 
#' phenotypes are also input for comparison and 
#' asseessment of predictor accuracy.
#'
#' @param inputdata output from filter_regions() \code{inputdata}
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
#' ## select regions to be used to build the predictor
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
	

test_predictor <- function(inputdata=NULL ,phenodata=NULL, phenotype=NULL, covariates=NULL,type="factor",predictordata=NULL){	
	requireNamespace("minfi", quietly=TRUE)
	requireNamespace("GenomicRanges", quietly=TRUE)
	requireNamespace("stats", quietly=TRUE)

	type <- match.arg(type,c("factor","binary", "numeric") )

	if(is.null(inputdata)) {
	 	stop('Must specify inputdata to use. This is the output from filter_regions()')
	}
	if(is.null(phenodata)) {
		stop('Must include phenotype file.')
	}
 	if(is.null(predictordata)) {
 		stop('Must specify predictor data to use. This is the output from build_predictor()')
 	}

 	predictor = predictordata
	
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
	regiondata = inputdata$regiondata[predictor$trainingProbes]


	ov = GenomicRanges::findOverlaps(inputdata$regiondata, predictor$regiondata)


	## predictions
	if(type=="factor"){
		# define possible predictions
		possibles = levels(droplevels(as.factor(phenodata[,phenotype])))
		possNA = c(possibles,"Unassigned")
		# make predictions
		minfi:::projectCellType(Y=expressiondata, coefCellType=predictor$coefEsts) -> predictions
		maxs <- apply(predictions,1,max)
		esttype = apply(predictions,1,which.highest)
		predicted <- possNA[esttype]  
	}
	if(type=="numeric"){
		## prepare data
		## ensure regions are named the same way as in build_predictor
		expressiondata = as.data.frame(t(expressiondata))
		colnames(expressiondata) <- paste0("exp_",1:ncol(expressiondata))

	
		knots_picked = predictor$knots_picked
		# Prepare model
		# Fit ns(expression, 5) for each expressed region
		l=5
		Xnew = model.matrix(as.formula(paste0("~",paste( paste0(" splines::ns(",colnames(expressiondata),",df=",l,", knots=knots_picked[,\'",colnames(knots_picked),"\'])"),collapse="+"))), data=expressiondata)

		## generate predictions
		predicted = as.numeric(as.matrix(t(predictor$coefEsts))%*% t(Xnew))
	}
	#summarize data
	actual <- phenodata[,phenotype]
	number_sites = length(predictor$trainingProbes)

	if(type=="factor"){
		number_match <- sum(predicted==actual)
		perc_correct = sum(predicted==actual)/length(actual)
		summarized = cbind(number_sites,number_match, perc_correct)
		colnames(summarized) <- c("sites_tested", "number_correct", "percent_correct")
	}
	if(type=="numeric"){
		correlation = stats::cor(predicted, actual)
		mean_diff = mean(abs(predicted-actual))
		summarized = cbind(number_sites, correlation, mean_diff)
		colnames(summarized) <- c("sites_tested", "correlation","mean_diff")
	}

	res <- list(actual = actual, predicted=predicted, summarized=summarized)

	return(res)

}