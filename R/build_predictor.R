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
#' @param type The class of the phenotype of interest (numeric, factor) 
#'\code{type}
#' @param covariates Which covariates to include in model \code{covariates}
#' @param numRegions The number of regions per class of variable of interest 
#' to pull out from each chromosome (default: 10) \code{numRegions}
#'
#' @return An n x m data.frame of coefficient estimates and region indices 
#' for each of the regions included from select_regions() along with regiondata and indices for trainingProbes
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
	require(caret)
	require(randomForest)
	require(gdata)

	## first, some checks
	 ## first, some checks
	type <- match.arg(type,c("factor", "numeric") )

	if(is.null(inputdata)) {
		stop('Must specify inputdata to use. This is the output from select_regions()')
	}
	if(is.null(phenodata)) {
		stop('Must include phenotype file.')
	}
	if(!(type %in% c('factor', 'binary', 'numeric'))) {
		stop('Phenotype you are predicting must be either "factor" or "numeric"')
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

	if(!is.null(covariates)){
	  	if(!all(covariates %in% names(pd))) {
	  		stop('Covariate included that is not in the prediction set. Please double check "covariates" argument.')
	 	 }
	  
		  ## pull out covariates to be included in the model
		  covar_data = pd[,covariates, drop=F]
		  ##drop unused levels for any factor covariates 
		  covar_data = gdata::drop.levels(covar_data)

		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	covars <- as.formula(paste("~ ", paste(covariates,collapse="+")))
	  	mm = model.matrix(covars, data=covar_data)
	  }
	
		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	message("Preparing Model")
	
	if(type=="factor"){  	
	  	## get list indeces for each group in the factor
		tIndexes <- split(seq_len(nrow(pd)), droplevels(pd[,phenotype, drop=F]))
		
				tstatList <- lapply(tIndexes, function(i) {
				    x <- rep(0, ncol(yGene))
				    x[i] <- 1 

				    if(!is.null(covariates)){
				  	  design = as.data.frame(cbind(x,mm))
				  	}else{
				  		design = as.data.frame(cbind(x = x))
				  	}  

			    fit = lmFit(yGene,design)			##### this is the SLOWWWW part
		    	eb = eBayes(fit)

		    	return(as.numeric(rownames(topTable(eb,1,n=numRegions))))
		    })
		      ## Note that in lmFit,
		      # g1mean <- rowMeans(normalized data in grp1)
		      # g2mean <- rowMeans(normalized data in grp2)
		      # fc <- g1mean - g2mean
			
		# in case not all have the number of probes
		cellSpecificList= lapply(tstatList, function(x) x[!is.na(x)])
		trainingProbes <- unique(unlist(cellSpecificList))
		regiondata = inputdata$regiondata[trainingProbes]

	}
	
	if(type=="numeric"){
		x=pd[,phenotype, drop=F]

	    if(!is.null(covariates)){
	  	  design = cbind(x = factor(x),mm)
	  	}else{
	  		design = cbind(x = factor(x))
	  	}

		fit = lmFit(yGene, design)
		eb = eBayes(fit)
		
		cellSpecificList = as.numeric(rownames(topTable(eb,1,n=numRegions)))
		trainingProbes = unique(cellSpecificList[!is.na(cellSpecificList)])
		regiondata = inputdata$regiondata[trainingProbes]
	}


	##################
	## Step2 
	# step 2: use the `minfi:::validationCellType` function to get the coefficients for prediction
	# you have to make the formula to pass, but we have some code in the `minfi:::pickCompProbes` function
	# it looks something like this, where `probeList` is really `cellSpecificInd` from above
	message("Calculating Coefficients")
	p=yGene
	
	## okay, now go ahead...
	p <- p[trainingProbes, ]
	pMeans <- colMeans(p)
	if(type=="factor"){
		pd[,phenotype] <- droplevels(as.factor(pd[,phenotype]))
		pd[,phenotype] <- gsub(" ","",pd[,phenotype])
		pd[,phenotype] <- factor(pd[,phenotype])
		names(pMeans) <- pd[,phenotype]
		form <- as.formula(sprintf("y ~ %s - 1", paste(levels(droplevels(as.factor(pd[,phenotype]))),
	    collapse = "+"))) 
		phenoDF <- as.data.frame(model.matrix(~pd[,phenotype] - 1))
		colnames(phenoDF) <- sub("^pd\\[, phenotype]", "", colnames(phenoDF))
		
		if (ncol(phenoDF) == 2) {
		    X <- as.matrix(phenoDF)
		    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
		}else{
		    tmp <- minfi:::validationCellType(Y = as.matrix(p), pheno = phenoDF, modelFix = form)
		    coefEsts <- tmp$coefEsts
		}

		res <- list(coefEsts = coefEsts, trainingProbes=trainingProbes, regiondata=regiondata)

	}
	if(type=="numeric"){

		# prepare data
		phenouse <- as.numeric(pd[,phenotype])
		datar = as.data.frame(t(p))
		colnames(datar) = names(regiondata)
		datar$phenouse = phenouse
		covars = pd[,covariates, drop=F]
		datar = cbind(datar, covars)
		
		# mdoel data
		plsFit <- train(phenouse ~ .,
                 data = datar,
                 method = "parRF",
                 ## Center and scale the predictors for the training
                ## set and all future samples.
                preProc = c("center", "scale"))
		res = list(coefEsts = plsFit, trainingProbes=trainingProbes, regiondata=regiondata)
	}

	return(res)
}