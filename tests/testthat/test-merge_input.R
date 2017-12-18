context("tests merge_input")

data(sysdata, package='phenopredict')
library(GenomicRanges)
## includes R object 'cm'   

## Make up some some region data
regions <- GRanges(seqnames = 'chr2', IRanges(
      start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
      end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))

## make up some expression data for 9 rows and 30 people
exp= cm[1:length(regions),1:30]

## Generate first object to be merged
inputdata1 <- list()
inputdata1$covmat = exp
inputdata1$regiondata = regions

## Generate scond object to be merged
regions2 = GRanges(seqnames = 'chr9', IRanges(
      start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
      end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))
exp2 = cm[1:length(regions2),1:30]
inputdata2 <- list()
inputdata2$covmat = exp2
inputdata2$regiondata = regions2

## merge objects
inputdata_merged<-phenopredict::merge_input(inputdata_list=list(inputdata1, inputdata2))

test_that("if merge_input output are correct class", {
  expect_is(inputdata_merged$regiondata, 'GRanges')
  expect_is(inputdata_merged$covmat, 'data.frame')
  expect_is(inputdata_merged, 'list')

})

test_that("output has expected dimensions", {
  expect_equal(length(inputdata_merged$regiondata), nrow(inputdata_merged$covmat))
})