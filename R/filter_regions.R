#' Filters regions from prediction set to be used to build the predctor for a 
#' phenotype of interest
#' 
#' This function takes phenotype of interest (sex, tissue type, etc.)
#' input by the user and uses a linear model (accounting for covariates, 
#' if provided) to filter those expressed regions that best predict the 
#' phenotype of interest. This is necessary when expression data are provided
#' in chunks or broken down by chromosome. These regions can then be merged together
#' with merge_regions() and are then used downstream for prediction. Default filters
#' top 100 expressed regions from input data.
#'
#' @param expression expression data where regions are in rows and samples are 
#' in columns \code{expression}
#' @param regiondata A GenomicRanges object in which  \code{regiondata}
#' @param phenodata phenotype data with samples in rows and corresponding 
#' phenotype 
#' information in columns  \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, factor)
#' \code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param numRegions The number of regions per class of variable of interest 
#' to pull out from each chromosome (default: 100) \code{numRegions}
#'
#' @return  The selected regions, the coverage matrix, and the region info 
#' to be used for prediction
#'
#' @keywords phenotype, prediction, filtering
#'
#' @importFrom stats model.matrix as.formula lm quantile
#' @importFrom gdata drop.levels  
#' @importFrom splines ns
#' @importFrom limma lmFit topTable eBayes
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
#' ## filter regions to be used to build the predictor
#' inputdata <- filter_regions(expression=exp, regiondata=regions,
#' 	phenodata=pheno, phenotype="sex", covariates=NULL,type="factor", numRegions=2)

filter_regions <- function(expression=NULL, regiondata=NULL ,phenodata=NULL, 
	phenotype=NULL, covariates=NULL,type=c("factor", "numeric"), numRegions=100){

	requireNamespace("limma", quietly=TRUE)
	requireNamespace("GenomicRanges", quietly=TRUE)
	requireNamespace("stats", quietly=TRUE)
	requireNamespace("gdata", quietly=TRUE)
	requireNamespace("splines", quietly=TRUE)

	# requireNamespace("genefilter", quietly=TRUE)

	## first, some checks
	type <- match.arg(type,c("factor", "numeric") )

	if(is.null(regiondata)) {
		stop('Must include a GenomicRanges object corresponding 
			to the regions included in expession')
	}
	if(is.null(expression)) {
		stop('Expression Data must be supplied.')
	}
	if(!(type %in% c('factor', 'numeric'))) {
		stop('Phenotype you are predicting must be either "factor" or "numeric"')
	}
	if(phenotype %in% covariates) {
		stop('Your phenotype of interest is also in your covariates. 
			Fix that first, please!')
	}
	 if(is.numeric(numRegions)==FALSE) {
	 	stop('Specify how many regions (per category type if 
	 		categorical) you want to filter with numRegions')
	  }

	  #### GET INDICES FOR PHENOTYPE OF INTEREST
	  yGene = expression
	  pd = phenodata
	  
	  if(!is.null(covariates)){
	  	if(!all(covariates %in% names(pd))) {
	  		stop('Covariate included that is not in the prediction set. 
	  			Please double check "covariates" argument.')
	 	 }
	  
		  ## pull out covariates to be included in the model
		  covar_data = as.data.frame(pd[,covariates, drop=FALSE])
		  covar_data = gdata::drop.levels(covar_data)

		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	covars <- stats::as.formula(paste("~ ", paste(covariates,collapse="+")))
	  	mm = stats::model.matrix(covars, data=covar_data)
	  }



	if(type=="factor"){  	
	  	## get list indeces for each group in the factor
		tIndexes <- split(seq_len(nrow(pd)), droplevels(pd[,phenotype, drop=F]))
		
				tstatList <- lapply(tIndexes, function(i) {
				    x <- rep(0, ncol(yGene))
				    x[i] <- 1 
				    				    x[tIndexes[[1]]] <- 1 

				    if(!is.null(covariates)){
				  	  design = as.data.frame(cbind(x,mm))
				  	}else{
				  		design = as.data.frame(cbind(x = x))
				  	}  
				## fit model for each level in phenotype
			    fit = limma::lmFit(yGene,design)			
			    ##### ^^this is the SLOWWWW part if you have a lot of regions
		    	eb = limma::eBayes(fit)

		    	return(as.numeric(rownames(limma::topTable(eb,coef=1,n=numRegions))))
		    })
		      ## Note that in lmFit,
		      # g1mean <- rowMeans(normalized data in grp1)
		      # g2mean <- rowMeans(normalized data in grp2)
		      # fc <- g1mean - g2mean
				
	cellSpecificList= lapply(tstatList, function(x) x[!is.na(x)])
	trainingProbes <- unique(unlist(cellSpecificList))

		# length(trainingProbes)
	}
	if(type=="numeric"){
		x=pd[,phenotype, drop=F]
		  
	    if(!is.null(covariates)){
	  		design = cbind(model.matrix(~splines::ns(get(phenotype),df=5)-1,data=pd),mm)
	  	}else{
	  		design = model.matrix(~splines::ns(get(phenotype),df=5)-1, data=pd )
	  	}

		fit = limma::lmFit(yGene,design)
		eb = limma::eBayes(fit)

		cellSpecificList = as.numeric(rownames(limma::topTable(eb,coef=1:5,n=numRegions)))
		trainingProbes = unique(cellSpecificList[!is.na(cellSpecificList)])

	}


		# just extract the regions we're going to use to build the predictor
		covmat = yGene[trainingProbes, ,drop=FALSE]
		regiondata = regiondata[trainingProbes, ,drop=FALSE]

		# make sure we know which sites those are
		index = c(as.numeric(trainingProbes))
		type=c(rep(c(phenotype),times=c(length(trainingProbes))))
		chr = GenomicRanges::seqnames(regiondata)	
		regioninfo = data.frame(chr=chr,type=type, 
		                  index=index)

		res <- list(regiondata=regiondata, covmat=covmat, regioninfo=regioninfo )
		return(res)



}
