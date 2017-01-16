#' Loads GTEx data to build predictor by chromosome  
#' 
#' This function (writtien by Leo Collado-Torres and Andrew Jaffee)
#' takes phenotype of interest (sex, tissue type, etc.)
#' input by the user and uses a linear model (accounting for covariates, 
#' if provided) to select those expressed regions that best predict the 
#' phenotype of interest. These regions are then used downstream for 
#' prediction. 
#'
#' @param chr chromosome to load (use UCSC naming scheme, ie: 'chr1') 
#' \code{chr}
#' @param minoverlap which output of annotateRegions() to load (1, 8, 20) 
#' \code{minoverlap}
#' @param chr_db  input database to use (ucsc, ensembl, gencode) \code{chr_db}
#' @param db database to use (ucsc, ensembl, gencode) \code{db}
#' @param phenoSize Specifies full or 'small' phenotype data to load 
#' \code{phenoSize}
#' @param help  specifies whether to surpress printing explanation message 
#' \code{help}
#' @param disk  either 'dcs' or 'dcl' \code{disk}
#'
#' @return res The GTEx data for the chromosome specified by input
#'
#' @keywords GTEx, expression, expressed region



gtexLoad <- function(chr = NULL, minoverlap = 1, chr_db = 'ucsc', db = 'ucsc', phenoSize = 'small', help = TRUE, disk = 'dcl') {
    require(GenomicRanges)
    require(GenomeInfoDb)

    ## added as.character twice before since which wasn't working using Jaffe's original function...when used in another function
    ## as.character(seqnames(regs_noEBV)
    stopifnot(db %in% c('ucsc', 'ensembl', 'gencode'))
    stopifnot(chr_db %in% c('ucsc', 'ensembl', 'gencode'))
    stopifnot(minoverlap %in% c(1, 8, 20))
    stopifnot(disk %in% c('dcl', 'dcs'))
    
    diskPath <- ifelse(disk == 'dcl', '/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis', '/dcs01/ajaffe/GTEX/Leek')
    
    regsFile <- paste0(diskPath, '/coverageMatrix/regions-cut0.5.Rdata')

    ## Select pheno file
    phenoFile <- paste0(diskPath, ifelse(phenoSize == 'small', '/pheno/pheno_missing_less_10.Rdata', '/pheno/pheno_complete.Rdata'))
    
    ## Read chr info
    gtexChr <- read.table(paste0(diskPath, '/coverageMatrix/simpleLoad/gtexChr.txt'), header = TRUE, colClasses = c('character', 'numeric', 'character', 'character', 'numeric', 'character', 'character', 'numeric'))
    
    if(disk == 'dcs') gtexChr$matrixFile <- gsub('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis', '/dcs01/ajaffe/GTEX/Leek', gtexChr$matrixFile)
    
    ## Select annotated regions file
    if(db == 'gencode') {
        annoFile <- NA
    } else {
        if(minoverlap != 1) {
            annoFile <- paste0(diskPath, '/coverageMatrix/annotatedRegions/annotated_', db, '_', minoverlap, '.Rdata')
        } else {
            annoFile <- paste0(diskPath, '/coverageMatrix/annotatedRegions/annotated_', db, '.Rdata')
        }
        
    }
    
    ## If chr is null, just return the info and the paths to the files
    if(is.null(chr)) {
        res <- list(chrInfo = gtexChr, regionsFile = regsFile, phenoFile = phenoFile, annotatedRegionsFile = annoFile)
        if(help) message(paste(Sys.time(), "Specify 'chr' (by default with the UCSC naming scheme) to load the data for that chromosome. The input naming scheme is determined by 'chr_db' and the output format by 'db'; each one has options 'ucsc', 'ensembl' and 'gencode'. 'phenoSize' == 'small' loads the pheno table where all variables have less than 10% missing, otherwise it loads the full pheno table. 'minoverlap' (1, 8 or 20) determines which output of annotateRegions() to load. 'disk' has to be either 'dcl' or 'dcs'. Using 'help' = FALSE disables printing this message when 'chr' is 'NULL'."))
        return(res)
    }
    
    ## Load regions and pheno data
    suppressPackageStartupMessages(library('GenomicRanges'))
    library('GenomeInfoDb')
    load(phenoFile)
    load(regsFile)

    ## Load data otherwise
    if(chr == 'chrEBV') {
        
        load(gtexChr$matrixFile[gtexChr$ucsc == 'chrEBV'])
        warning("'db' is ignored since chrEBV is not part of the UCSC knownGene or Ensembl annotations. There won't be any annotated regions.")
        res <- list(chrInfo = subset(gtexChr, ucsc == 'chrEBV'), coverageMatrix = coverageMatrix, regions = regions[seqnames(regions) == 'chrEBV'], pheno = pheno, annotatedRegions = NULL)
        return(res)
    }

    
    stopifnot(chr %in% gtexChr[, chr_db])

    ## Subset regions
    chr_ucsc <- gtexChr$ucsc[gtexChr[, chr_db] == chr]
    chr_ucsc <- chr_ucsc[!is.na(chr_ucsc)]
    regs <- regions[seqnames(regions) == chr_ucsc]
    regs <- keepSeqlevels(regs, chr_ucsc)
    seqlevels(regs) <- gtexChr[gtexChr$ucsc == chr_ucsc, db]

    ## Load annotated regions info
    if(is.null(annoFile)) {
        annoRegs <- NULL
    } else {
        load(annoFile)
        regs_noEBV <- dropSeqlevels(regions, 'chrEBV')

        assign('annoRegions', get(ls()[grep('annotated_', ls())]))
        gs_present <- names(annoRegions$annotationList)
    
        annoRegs <- list(
            countTable = annoRegions$countTable[as.vector(seqnames(regs_noEBV) == chr_ucsc), ],
            annotationList = annoRegions$annotationList[gs_present %in% as.character(which(as.character(seqnames(regs_noEBV)) == chr_ucsc))]
        )
        names(annoRegs$annotationList) <- as.integer(names(annoRegs$annotationList)) - (which(as.character(seqnames(regs_noEBV)) == chr_ucsc)[1] - 1)
    }
    
    ## Load coverage matrix if it's available
    matFile <- gtexChr$matrixFile[gtexChr$ucsc == chr_ucsc]
    if(!is.na(matFile)) {
        load(matFile)
    } else {
        coverageMatrix <- NULL
    }
    
    res <- list(chrInfo = subset(gtexChr, ucsc == chr_ucsc), coverageMatrix = coverageMatrix, regions = regs, pheno = pheno, annotatedRegions = annoRegs)
    return(res)
        
}
