#' Selects regions from prediction set to be used to build the predctor for a 
#' phenotype of interest
#' 
#' This function takes phenotype of interest (sex, tissue type, etc.)
#' input by the user and uses a linear model (accounting for covariates, 
#' if provided) to select those expressed regions that best predict the 
#' phenotype of interest. These regions are then used downstream for 
#' prediction. 
#'
#' @param expression expression data where regions are in rows and samples are 
#' in columns \code{expression}
#' @param regiondata A GRanges object in which  \code{regiondata}
#' @param phenodata phenotype data with samples in rows and corresponding 
#' phenotype 
#' information in columns  \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, binary, factor)
#' \code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param numRegions The number of regions per class of variable of interest 
#' to pull out from each chromosome (default: 100) \code{numRegions}
#'
#' @return  The selected regions, the coverage matrix, and the region info 
#' to be used for prediction
#'
#' @keywords phenotype, prediction, selection
#'
#' @export
#' 

select_regions <- function(expression=NULL, regiondata=NULL ,phenodata=NULL, phenotype=NULL, covariates=NULL,type=c("factor","binary","numeric"), numRegions=100){

	require(limma)
	require(GenomicRanges)
	require(stats)

	## first, some checks
	type <- match.arg(type)

	if(is.null(regiondata)) {
		stop('Must include a GRanges object corresponding to the regions included in expession')
	}
	if(is.null(expression)) {
		stop('Expression Data must be supplied.')
	}
	if(!(type %in% c('factor', 'binary', 'factor'))) {
		stop('Phenotype you are predicting must be either "factor","binary", or "numeric"')
	}
	if(phenotype %in% covariates) {
		stop('Your phenotype of interest is also in your covariates. Fix that first, please!')
	}
	if(is.numeric(numRegions)==FALSE) {
		stop('Specify how many regions per category type you want to select with numRegions')
	}

	  #### GET INDICES FOR PHENOTYPE OF INTEREST
	  yGene = expression
	  pd = phenodata
	  if(!all(covariates %in% names(pd))) {
	  	stop('Covariate included that is not in the prediction set. Please double check "covariates" argument.')
	  }
	  
	  ## pull out covariates to be included in the model
	  covars = pd[,covariates]
	  ##get covariates in order
	

		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	covars <- as.formula(paste("~ ", paste(covariates,collapse="+")))
	  	mm = model.matrix(covars, data=pd)

	  	## get list indeces for each group in the factor
		tIndexes <- split(seq_len(pd), droplevels(pd[,phenotype]))
		tstatList <- lapply(tIndexes, function(i) {
		    x <- rep(0, ncol(yGene))
		    x[i] <- 1 		

		    design = cbind(x = x,mm)
		     
		    fit = lmFit(yGene,design)			##### this is the SLOWWWW part
		    eb = eBayes(fit)
		    return(as.numeric(rownames(topTable(eb,1,n=numRegions))))
		      ## Note that in lmFit,
		      # g1mean <- rowMeans(normalized data in grp1)
		      # g2mean <- rowMeans(normalized data in grp2)
		      # fc <- g1mean - g2mean

		})

		# in case not all have the number of probes
		cellSpecificList= lapply(tstatList, function(x) x[!is.na(x)])
		trainingProbes <- unique(unlist(cellSpecificList))

		# just extract the regions we're going to use to build the predictor
		covmat = yGene[trainingProbes,]
		regiondata = regiondata[trainingProbes]

		# make sure we know which sites those are
		index = c(as.numeric(trainingProbes))
		type=c(rep(c(phenotype),times=c(length(trainingProbes))))
		chr = seqnames(regiondata)	
		regioninfo = data.frame(chr=chr,type=type, 
		                  index=index)

		res <- list(regiondata=regiondata, covmat=covmat, regioninfo=regioninfo )
		return(res)


}