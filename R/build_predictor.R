#' Uses output from select_regions() to build the predictor   
#' 
#' This function takes the output from select_regions() and builds the 
#' phenotype predictor to be used for phenotype prediction in a new 
#' data set. This function outputs the fits (coefficient estimates, `coefEsts`)
#' from the model for the selectied regions as well as the rows selected 
#' (`trainingProbes`) from the input data.  
#'
#' @param inputdata output from select_regions() \code{inputdata}
#' @param phenodata data set with phenotype information; samples in rows, 
#' variables in columns \code{phenodata}
#' @param phenotype phenotype of interest \code{phenotype}
#' @param type The class of the phenotype of interest (numeric, binary, factor 
#'\code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param numRegions The number of regions per class of variable of interest 
#' to pull out from each chromosome (default: 10) \code{numRegions}
#'
#' @return An n x m data.frame of coefficient estimates and region indices 
#' for each of the regions included from select_regions()
#'
#' @keywords phenotype, prediction, coefficient estimates
#'
#' @export
#' 


build_predictor <- function(inputdata=NULL ,phenodata=NULL, phenotype=NULL, covariates=NULL,type=NULL, numRegions=10){
	require(GenomicRanges)
	require(limma)
	require(stats)
	require(minfi)

	## first, some checks
	 if(is.null(type)) {
	  	stop('Must specify which type of phenotype you are interested in predicting ("factor","binary","numeric")')
	  }
	  if(is.null(inputdata)) {
	  	stop('Must specify inputdata to use. This is the output from select_regions()')
	  }
 	 if(is.null(phenodata)) {
	  	stop('Must include phenotype file.')
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
	  if(ncol(inputdata$covmat) != nrow(phenodata)) {
	  	stop('The number of samples in your inputdata must be the same as in your phenotype file')
	  }



	##################
	  yGene = inputdata$covmat
	  pd = phenodata

	  if(!all(covariates %in% names(pd))) {
	  	stop('Covariate included that is not in the prediction set. Please double check "covariates" argument.')
	  	# stopifnot(covariates %in% names(pd))
	  }
	
		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	message("Preparing Model")
	  	covars <- as.formula(paste("~ ", paste(covariates,collapse="+")))
	  	mm = model.matrix(covars, data=pd)
	  	
		tIndexes <- split(1:nrow(pd), droplevels(pd[,phenotype]))
		tstatList <- lapply(tIndexes, function(i) {
		    x <- rep(0, ncol(yGene))
		    x[i] <- 1

		    # design = model.matrix(~factor(x) + covars,data=pd)
		    design = cbind(x = x,mm)
		      fit = lmFit(yGene,design)
		      eb = eBayes(fit)
		      return(topTable(eb,1,n=numRegions))
		 })

	cellSpecificList <- lapply(tstatList, function(x) {
	    y <- x[x[,"P.Value"]< 1e-05,]
	    values = abs(y[,"logFC"]) ## code changed to allow for tissue expression to be higher or lower than expression in other tissues
	    y <- y[order(values,decreasing=TRUE),] # negatives are higher in cell type
	    c(rownames(y)[1:numRegions])
	})
	# in case not all have the number of probes
	cellSpecificList= lapply(cellSpecificList, function(x) x[!is.na(x)])
	cellSpecificInd = lapply(cellSpecificList,
	    function(x) rownames(yGene) %in% x)



	##################
	## Step2 
	# step 2: use the `minfi:::validationCellType` function to get the coefficients for prediction
	# you have to make the formula to pass, but we have some code in the `minfi:::pickCompProbes` function
	# it looks something like this, where `probeList` is really `cellSpecificInd` from above
	message("Calculating Coefficients")
	p=yGene
	probeList=cellSpecificList
	trainingProbes <- as.numeric(unique(unlist(probeList)))

	## okay, now go ahead...
	p <- p[trainingProbes, ]
	pMeans <- colMeans(p)
	pd[,phenotype] <- droplevels(pd[,phenotype])
	pd[,phenotype] <- gsub(" ","",pd[,phenotype])
	pd[,phenotype] <- factor(pd[,phenotype])
	names(pMeans) <- pd[,phenotype]
	form <- as.formula(sprintf("y ~ %s - 1", paste(levels(droplevels(pd[,phenotype])),
	    collapse = "+"))) 
	phenoDF <- as.data.frame(model.matrix(~pd[,phenotype] - 1))
	colnames(phenoDF) <- sub("^pd[,phenotype]", "", colnames(phenoDF))
	if (ncol(phenoDF) == 2) {
	    X <- as.matrix(phenoDF)
	    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
	}else{
	    tmp <- minfi:::validationCellType(Y = p, pheno = phenoDF, modelFix = form)
	    coefEsts <- tmp$coefEsts
	}

	res <- list(coefEsts = coefEsts, trainingProbes=trainingProbes)

	return(res)
}