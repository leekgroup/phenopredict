#' Example epression data used to build the predictor.
#'
#' A dataset containing expression data for 500 expression regions
#' on chrX (generated from GTEx data)
#' 
#' @format A data frame with 500 rows and 4769 samples (columns):
#'
#' @source \url{https://bioconductor.org/packages/release/bioc/html/recount.html}
"cm"

#' Example expression data to be used for prediction.
#'
#' A dataset containing expression data for 500 expression regions
#' on chrX (generated from GTEx data)
#' 
#' @format A data frame with 500 rows and 4769 different samples (columns)
#' than those used in "cm":
#'
#' @source \url{https://bioconductor.org/packages/release/bioc/html/recount.html}
"cm_new"

#' Example phenotype data.
#'
#' A dataset containing phenotype data for 4769 samples included in "cm"
#' 
#' @format A data frame with 4769 rows and 109 columns
#'
#' @source \url{https://bioconductor.org/packages/release/bioc/html/recount.html}
"pheno"

#' Information about regions included
#'
#' A GRanges object with information for the 500 expressed regions included in cm
#' and cm_new
#' 
#' @format A GRanges object with 500 ranges
"regiondata"

#' Output from select_regions
#'
#' A list including a GRanges object (regiondadta), a expression matrix (cm),
#' and a dataframe (regioninfo). This is the output of select_regions().
#' 
#' @format A list containing three elements
#' \describe{
#'   \item{regioninfo}{data.matrix, region information}
#'   \item{regiondata}{GRanges, genomic ranges region information}
#'   \item{covmat}{matrix, expression coverage}
#' }
"inputdata"

#' Output from build_predictor
#'
#' A list including a data.frame of coefficient estimates
#' and a vector of probes used in the model. This is the 
#' output of build_predictor().
#' 
#' @format A list containing two elements
#' \describe{
#'   \item{coefEsts}{data.frame, prediction estimates}
#'   \item{trainingProbes}{vector, probes to be extracted}
#' }
"predictor"