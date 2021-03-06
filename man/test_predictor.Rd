% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_predictor.R
\name{test_predictor}
\alias{test_predictor}
\title{Test accuracy of predictor on known phenotypes}
\usage{
test_predictor(inputdata = NULL, phenodata = NULL, phenotype = NULL,
  covariates = NULL, type = "factor", predictordata = NULL)
}
\arguments{
\item{inputdata}{output from filter_regions() \code{inputdata}}

\item{phenodata}{data set with phenotype information; samples in rows,
variables in columns \code{phenodata}}

\item{phenotype}{phenotype of interest \code{phenotype}}

\item{covariates}{Which covariates to include in model \code{covariates}}

\item{type}{The class of the phenotype of interest (numeric, binary, factor)
\code{type}}

\item{predictordata}{object output from build_predictor \code{predictordata}}
}
\value{
list of actual and predicted phenotype, and summarization of output
}
\description{
This function takes the expression data input to
build_predictor() and the coefficient estimates from
build_predictor() for phenotype prediction. The known
phenotypes are also input for comparison and
asseessment of predictor accuracy.
}
\examples{

library('GenomicRanges')
library('dplyr')

## Make up some some region data
regions <- GRanges(seqnames = 'chr2', IRanges(
     start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
     end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))

## make up some expression data for 9 rows and 30 people
data(sysdata, package='phenopredict')
## includes R object 'cm'
exp= cm[1:length(regions),1:30]

## generate some phenotype information
sex = as.data.frame(rep(c("male","female"),each=15))
age = as.data.frame(sample(1:100,30))
pheno = dplyr::bind_cols(sex,age)
colnames(pheno) <- c("sex","age")

## select regions to be used to build the predictor
inputdata <- filter_regions(expression=exp, regiondata=regions,
	phenodata=pheno, phenotype="sex",
	covariates=NULL,type="factor", numRegions=2)

## build phenotype predictor
predictor<-build_predictor(inputdata=inputdata ,phenodata=pheno,
	phenotype="sex", covariates=NULL,type="factor", numRegions=2)

## determine resubstitution error
## carry out prediction in training data set
predictions_test<-test_predictor(inputdata=inputdata ,phenodata=pheno,
	phenotype="sex", covariates=NULL,type="factor",predictordata=predictor)
}
\keyword{phenotype,}
\keyword{prediction,}
\keyword{test}
