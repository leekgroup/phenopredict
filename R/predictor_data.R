#' Extracts expression info from prediction set (GTEx) using output from select_regions()  
#' 
#' This function takes uses the regions selected from select_regions and extracts the
#' expression and region information necessary to build the predictor
#'
#' @param phenotype phenotype of interest \code{phenotype}
#' @param chrom Which chromosomes to include in predictor \code{chrom}
#'
#' @return output Three files: regions_to_subset, coverageMatrixPredictor, and regioninfo
#'
#' @keywords phenotype, prediction, selection
#'
#' @export
#' 
#' @examples
#' predictor_data(phenotype="SMTS",chrom=22)

predictor_data <- function(phenotype=NULL,chrom="ALL"){
	
	##########
	### Get the unique regions
	### from the GTEX selection proces
	##########
	
	### Load GTEX loading function
	source("/users/sellis/phenopredict/R/gtexLoad.R")

	### Loop over chromosomes and load data, fitting the models
	if(chrom=="ALL"){
		chrs = c("X","Y",1:22)
	}else{
		chrs=chrom
	}

	regions_to_subset = GRanges()
	for(i in 1:length(chrs)){
	  chrname=paste0("chr",chrs[i])
	  dat = gtexLoad(chrname)
	  load(paste0(chrname,"_",phenotype,".rda"))
	  regions_to_subset <- append(regions_to_subset, dat$regions[unique(out$index)])
	  rm(out)
	  cat(i)
	}
	save(regions_to_subset,file=paste0("regions_to_subset_",phenotype,".rda"))

	#####
	### Save GTEX coverage
	#####
	## only includes the samples used to build the predictor
	covmat = matrix(NA,nrow=1,ncol=ncol(dat$coverageMatrix))
	for(i in 1:length(chrs)){
	  chrname=paste0("chr",chrs[i])
	  dat = gtexLoad(chrname)
	  load(paste0(chrname,"_",phenotype,".rda"))
	  covmat = rbind(covmat,dat$coverageMatrix[unique(out$index),])
	  cat(i)
	}
	covmat = covmat[-1,]
	covmat <- covmat[,usegtex]
	covmat <- covmat[,sample_individuals]
	coverageMatrixPredictor = covmat
	save(coverageMatrixPredictor,file=paste0("coverageMatrixPredictor_",phenotype,".rda"))

	#### 
	## Save region information
	####	
	regioninfo = data.frame(type=NA,index=NA,chr=NA)
	for(i in 1:length(chrs)){
	  chrname=paste0("chr",chrs[i])
	  load(paste0(chrname,"_",phenotype,".rda"))
	  out$chr = chrname
	  regioninfo = rbind(regioninfo,out)
	  cat(i)
	}
	regioninfo=regioninfo[-1,]
	save(regioninfo,file=paste0("regioninfo_",phenotype,".rda"))
}