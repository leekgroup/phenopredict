#' Uses output from filter_regions() to build the predictor   
#' 
#' This function takes the output from filter_regions() and builds the 
#' phenotype predictor to be used for phenotype prediction in a new 
#' data set. This function outputs the fits (coefficient estimates, `coefEsts`)
#' from the model for the selectied regions as well as the rows selected 
#' (`trainingProbes`) from the input data.  
#'
#' @param inputdata output from filter_regions() \code{inputdata}
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
#' for each of the regions included from filter_regions() 
#' along with regiondata and indices for trainingProbes
#'
#' @keywords phenotype, prediction, coefficient estimates
#'
#' @importFrom stats model.matrix as.formula lm quantile
#' @importFrom gdata drop.levels  
#' @importFrom splines ns
#' @importFrom limma lmFit topTable eBayes
#' @importFrom broom tidy
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
#' phenodata=pheno, phenotype="sex", covariates=NULL,type="factor", 
#' numRegions=2)
#' 
#' ## build phenotype predictor
#' predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno, 
#' 	phenotype="sex", covariates=NULL,type="factor", numRegions=2)

build_predictor <- function(inputdata=NULL ,phenodata=NULL, 
	phenotype=NULL, covariates=NULL,type=NULL, numRegions=NULL){
	requireNamespace("GenomicRanges", quietly=TRUE)
	requireNamespace("limma", quietly=TRUE)
	requireNamespace("stats", quietly=TRUE)
	requireNamespace("minfi", quietly=TRUE)
	requireNamespace("gdata", quietly=TRUE)
	requireNamespace("splines", quietly=TRUE)
	requireNamespace("broom", quietly=TRUE)

	## first, some checks
	type <- match.arg(type,c("factor", "numeric") )

	if(is.null(inputdata)) {
		stop('Must specify inputdata to use. 
			This is the output from filter_regions()')
	}
	if(is.null(phenodata)) {
		stop('Must include phenotype file.')
	}
	if(!(type %in% c('factor', 'binary', 'numeric'))) {
		stop('Phenotype you are predicting must be either "factor" or "numeric"')
	}
	if(phenotype %in% covariates) {
		stop('Your phenotype of interest is also in your covariates. 
			Fix that first, please!')
	}
	if(is.numeric(numRegions)==FALSE) {
	 	stop('Specify how many regions per category type you want 
	 		to select with numRegions')
	 }
	if(ncol(inputdata$covmat) != nrow(phenodata)) {
		stop('The number of samples in your inputdata must be the 
			same as in your phenotype file')
	}



	##################
	yGene = inputdata$covmat
	pd = phenodata

	if(!is.null(covariates)){
	  	if(!all(covariates %in% names(pd))) {
	  		stop('Covariate included that is not in the prediction set. 
	  			Please double check "covariates" argument.')
	 	 }
	  
		  ## pull out covariates to be included in the model
		  covar_data = pd[,covariates, drop=FALSE]
		  ##drop unused levels for any factor covariates 
		  covar_data = gdata::drop.levels(covar_data)

		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	covars <- stats::as.formula(paste("~ ", paste(covariates,collapse="+")))
	  	mm = stats::model.matrix(covars, data=covar_data)
	  }
	
		## instead of using rowttests, use LmFit 
		## to compute differences & to include covariates
	  	# message("Preparing Model")
	
	if(type=="factor"){  	
	  	## get list indeces for each group in the factor
		tIndexes <- split(seq_len(nrow(pd)), droplevels(pd[,phenotype, drop=FALSE]))
		
				tstatList <- lapply(tIndexes, function(i) {
				    x <- rep(0, ncol(yGene))
				    x[i] <- 1 

				    if(!is.null(covariates)){
				  	  design = as.data.frame(cbind(x,mm))
				  	}else{
				  		design = as.data.frame(cbind(x = x))
				  	}  

			    fit = limma::lmFit(yGene,design)			##### this is the SLOWWWW part
		    	eb = limma::eBayes(fit)

		    	return(as.numeric(rownames(limma::topTable(eb,1,n=numRegions))))
		    })
		      ## Note that in lmFit,
		      # g1mean <- rowMeans(normalized data in grp1)
		      # g2mean <- rowMeans(normalized data in grp2)
		      # fc <- g1mean - g2mean
			
		# in case not all have the number of probes
		cellSpecificList= lapply(tstatList, function(x) x[!is.na(x)])
		trainingProbes <- unique(unlist(cellSpecificList))
	}
	
	if(type=="numeric"){
		x=pd[,phenotype, drop=FALSE]
		  
	    if(!is.null(covariates)){
	  		design = cbind(model.matrix(~splines::ns(
	  			get(phenotype),df=5) - 1,data=pd),mm)
	  	}else{
	  		design = model.matrix(~splines::ns(get(phenotype),df=5) - 1, data=pd)
	  	}

		fit = limma::lmFit(yGene,design)
		eb = limma::eBayes(fit)

		cellSpecificList = as.numeric(rownames(limma::topTable(eb,1:5,n=numRegions)))
		trainingProbes = unique(cellSpecificList[!is.na(cellSpecificList)])
	}

	regiondata = inputdata$regiondata[trainingProbes]



	##################
	## Step2 
	# step 2: use the `minfi::validationCellType` function 
	# to get the coefficients for prediction
	# you have to make the formula to pass, but we have some code 
	# in the `minfi:::pickCompProbes` function
	# it looks something like this, where `probeList` is really 
	# `cellSpecificInd` from above
	# message("Calculating Coefficients")
	p=yGene
	
	## okay, now go ahead...
	p <- p[trainingProbes, ]
	# pMeans <- colMeans(p)
	if(type=="factor"){
		pd[,phenotype] <- droplevels(as.factor(pd[,phenotype]))
		pd[,phenotype] <- gsub(" ","",pd[,phenotype])
		pd[,phenotype] <- factor(pd[,phenotype])
		# names(pMeans) <- pd[,phenotype]
		form <- stats::as.formula(sprintf("y ~ %s - 1", 
			paste(levels(droplevels(as.factor(pd[,phenotype]))),
	    collapse = "+"))) 
		phenoDF <- as.data.frame(stats::model.matrix(~pd[,phenotype] - 1))
		colnames(phenoDF) <- sub("^pd\\[, phenotype]", "", colnames(phenoDF))
		
		if (ncol(phenoDF) == 2) {
		    X <- as.matrix(phenoDF)
		    coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
		}else{
		    tmp <- minfi:::validationCellType(Y = as.matrix(p), 
		    	pheno = phenoDF, modelFix = form)
		    coefEsts <- tmp$coefEsts
		}

		res <- list(coefEsts = coefEsts, trainingProbes=trainingProbes, 
			regiondata=regiondata)

	}
	if(type=="numeric"){
		
		# in build predictor do this
		## check to make sure that 5*ER !> N throw an error, 
		## you have more expressed region terms w/ spline than sample size
		exp_data = as.data.frame(t(p[trainingProbes, ]))
		colnames(exp_data) <- paste0("exp_",1:ncol(exp_data))

		# Prepare model
		# Fit ns(expression, 5) for each expressed region
		l=5

		## set knots to values that are nonzero
		knot_pick <- function(x){
			## determine quantile values for expression
			
			vals <- seq(0,1,0.1)
			quants <- quantile(x, probs=vals)
			#figure out where to put first knot
			minval <- min(which(quants > 0))
			# space evenly from there 
			if(minval==1){
				knot_val <- c(0.2,0.4,0.6,0.8)
				#Then the basis will include two boundary knots and 4 internal knots, 
				# placed at the 20th, 40th, 60th, and 80th quantiles of x, respectively. 
				# The boundary knots, by default, are placed at the min and max of x.

			}else{
				spacing <- (1-vals[minval])/4
				knot_val <- seq(vals[minval],1,spacing)[1:4]
			}
			# define knots
			knots_x <- quantile(x, probs=knot_val)
		}	
		knots_picked <- apply(exp_data,2,knot_pick)


		
		 X = model.matrix(as.formula(paste0("~",
		 	paste( paste0(" splines::ns(",colnames(exp_data),",df=",l,", 
		 	knots=knots_picked[,\'",colnames(knots_picked),"\'])"),
		 	collapse="+"))), data=exp_data)

		## Model data
		lm1 = lm(pd[,phenotype] ~ X,data=exp_data)

		## get coefEsts
		coefEsts = broom::tidy(lm1)[,2]

		#get coefficient estimates for regions included in topTable		
		res <- list(coefEsts = coefEsts, trainingProbes=trainingProbes, 
			regiondata=regiondata, knots_picked=knots_picked)

	}

	return(res)
}