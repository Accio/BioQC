## test wmwTest
library(BioQC)

context("Non-exact wilcoxon test in R")
testNums <- 1:10
testSub <- rep_len(c(TRUE, FALSE), length.out=length(testNums))
testP <- wmwTestInR(testNums, testSub)
testU <- wmwTestInR(testNums, testSub, valType="U")
testIndP <- wmwTestInR(testNums, which(testSub))
testIndU <- wmwTestInR(testNums, which(testSub), valType="U")
expWMW <- wilcox.test(testNums[seq(1,9,2)],
                      testNums[seq(2,10,2)], exact=FALSE)

test_that("Non-exact wilcoxon test in R is wrapped correctly by wmwTestInR", {
             expect_equal(testP, expWMW$p.value)
             expect_equal(testIndP, expWMW$p.value)
             expect_equivalent(testU, expWMW$statistic)
             expect_equivalent(testIndU, expWMW$statistic)
         })

context("Wilcoxon test in R and C")

set.seed(1887)
NROW <- 2000
NCOL <- 5
GSCOUNT <- 50
TOL <- 1E-2
GSSIZE <- sample(1:as.integer(NROW/2), GSCOUNT, replace=TRUE)
ind <- lapply(GSSIZE, function(x) sample(1:NROW, x))
exprs <- matrix(round(rnorm(NROW*NCOL),4), nrow=NROW)

gc()
system.time(Cless <- wmwTest(exprs, ind, valType="p.less"))
system.time(Cgreater <- wmwTest(exprs, ind, valType="p.greater"))
system.time(Ctwosided <- wmwTest(exprs, ind, valType="p.two.sided"))
wmwTestR <- function(matrix, index, valType) {
  sub <- rep(FALSE, length(matrix))
  sub[index]=TRUE
  wmwTestInR(matrix, sub, valType)
}
system.time(Rless <- apply(exprs, 2, function(x)
                          sapply(ind, function(y) 
                                 wmwTestR(x, y, valType="p.less"))))
system.time(Rgreater <- apply(exprs, 2, function(x)
                              sapply(ind, function(y)
                                     wmwTestR(x, y, valType="p.greater"))))
system.time(Rtwosided <- apply(exprs, 2, function(x)
                              sapply(ind, function(y)
                                     wmwTestR(x, y, valType="p.two.sided"))))

test_that("less is identical", {
              expect_equal(Cless, Rless, tolerance=TOL)
          })
test_that("greater is identical", {
              expect_equal(Cgreater, Rgreater,tolerance=TOL)
          })
test_that("two.sided is identical", {
              isNearOne <- 1-Ctwosided<0.01
              expect_equal(object=Ctwosided[!isNearOne], expected=Rtwosided[!isNearOne], tolerance=TOL, scale=1L)
          })


## TO BE REFACTORED!
context("one-to-one copy of the example")
## R-native data structures
set.seed(1887)
rd <- rnorm(1000)
rl <- sample(c(TRUE, FALSE), 1000, replace=TRUE)
wmwTest(rd, rl, valType="p.two.sided")
wmwTest(rd, which(rl), valType="p.two.sided")
rd1 <- rd + ifelse(rl, 0.5, 0)
wmwTest(rd1, rl, valType="p.greater")
wmwTest(rd1, rl, valType="U")
rd2 <- rd - ifelse(rl, 0.2, 0)
wmwTest(rd2, rl, valType="p.greater")
wmwTest(rd2, rl, valType="p.two.sided")
wmwTest(rd2, rl, valType="p.less")

## matrix forms
rmat <- matrix(c(rd, rd1, rd2), ncol=3, byrow=FALSE)
wmwTest(rmat, rl, valType="p.two.sided")
wmwTest(rmat, rl, valType="p.greater")

wmwTest(rmat, which(rl), valType="p.two.sided")
wmwTest(rmat, which(rl), valType="p.greater")

## other valTypes
wmwTest(rmat, which(rl), valType="U")
wmwTest(rmat, which(rl), valType="abs.log10p.greater")
wmwTest(rmat, which(rl), valType="log10p.less")
wmwTest(rmat, which(rl), valType="abs.log10p.two.sided")
wmwTest(rmat, which(rl), valType="Q")

## using ExpressionSet
data(sample.ExpressionSet)
testSet <- sample.ExpressionSet
fData(testSet)$GeneSymbol <- paste("GENE_",1:nrow(testSet), sep="")
mySig1 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.25, 0.75), replace=TRUE)
wmwTest(testSet, which(mySig1), valType="p.greater")

## using integer
exprs(testSet)[,1L] <- exprs(testSet)[,1L] + ifelse(mySig1, 50, 0)
wmwTest(testSet, which(mySig1), valType="p.greater")

## using lists
mySig2 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.6, 0.4), replace=TRUE)
wmwTest(testSet, list(first=mySig1, second=mySig2))

## using GMT file
gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt_list <- readGmt(gmt_file)

gss <- sample(unlist(sapply(gmt_list, function(x) x$genes)), 1000)
eset<-new("ExpressionSet",
          exprs=matrix(rnorm(10000), nrow=1000L),
          phenoData=new("AnnotatedDataFrame", data.frame(Sample=LETTERS[1:10])),
          featureData=new("AnnotatedDataFrame",data.frame(GeneSymbol=gss)))
esetWmwRes <- wmwTest(eset ,gmt_list, valType="p.greater")
summary(esetWmwRes)
        
