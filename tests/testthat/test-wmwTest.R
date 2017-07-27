## test wmwTest
library(BioQC)

context("Non-exact wilcoxon test in R")
testNums <- 1:10
testSub <- rep_len(c(TRUE, FALSE), length.out=length(testNums))
testP <- wmwTestInR(testNums, testSub, valType="p.two.sided")
testW <- wmwTestInR(testNums, testSub, valType="W")
testIndP <- wmwTestInR(testNums, which(testSub), valType="p.two.sided")
testIndW <- wmwTestInR(testNums, which(testSub), valType="W")
expWMW <- wilcox.test(testNums[seq(1,9,2)],
                      testNums[seq(2,10,2)], alternative="two.sided", exact=FALSE)

test_that("Non-exact wilcoxon test in R is wrapped correctly by wmwTestInR", {
             expect_equal(testP, expWMW$p.value)
             expect_equal(testIndP, expWMW$p.value)
             expect_equivalent(testW, expWMW$statistic)
             expect_equivalent(testIndW, expWMW$statistic)
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


context("Test different combinations of input and valType in wmwTest")
## R-native data structures

rd <- rnorm(1000)
rl <- sample(c(TRUE, FALSE), 1000, replace=TRUE)
wmwTwoSidedC <- wmwTest(rd, rl, valType="p.two.sided")
wmwTwoSidedR <- wmwTestInR(rd, rl, valType="p.two.sided")

rd1 <- rd + ifelse(rl, 0.5, 0)
wmwGreaterC <- wmwTest(rd1, rl, valType="p.greater")
wmwGreaterR <- wmwTestInR(rd1, rl, valType="p.greater")

rd2 <- rd - ifelse(rl, 0.2, 0)
wmwGreaterC2 <- wmwTest(rd2, rl, valType="p.greater")
wmwTwoSidedC2 <- wmwTest(rd2, rl, valType="p.two.sided")
wmwLessC2 <- wmwTest(rd2, rl, valType="p.less")
wmwGreaterR2 <- wmwTestInR(rd2, rl, valType="p.greater")
wmwTwoSidedR2 <- wmwTestInR(rd2, rl, valType="p.two.sided")
wmwLessR2 <- wmwTestInR(rd2, rl, valType="p.less")
test_that("p-values by C and R implementations are identical", {
    expect_equivalent(wmwTwoSidedC, wmwTwoSidedR)
    expect_equivalent(wmwGreaterC, wmwGreaterR)
    expect_equivalent(wmwGreaterC2, wmwGreaterR2)
    expect_equivalent(wmwTwoSidedC2, wmwTwoSidedR2)
    expect_equivalent(wmwLessC2, wmwLessR2)
})

wmwTwoSidedC.ind <- wmwTest(rd, which(rl), valType="p.two.sided")
test_that("wmwTest accepts logical and numeric vector as IndexList", {
    expect_equivalent(wmwTwoSidedC.ind, wmwTwoSidedC)
})

## matrix forms
rmat <- matrix(c(rd, rd1, rd2), ncol=3, byrow=FALSE)
tsMat <- wmwTest(rmat, rl, valType="p.two.sided")
tsMatR <- apply(rmat, 2, function(x)
                wmwTestInR(x, rl, valType="p.two.sided"))
tsMatLess <- wmwTest(rmat, rl, valType="p.less")
tsMatLessR <- apply(rmat, 2, function(x)
                       wmwTestInR(x, rl, valType="p.less"))
tsMatGreater <- wmwTest(rmat, rl, valType="p.greater")
tsMatGreaterR <- apply(rmat, 2, function(x)
                       wmwTestInR(x, rl, valType="p.greater"))
logicTsMat <- wmwTest(rmat, which(rl), valType="p.two.sided")
logicTsMatGreater <- wmwTest(rmat, which(rl), valType="p.greater")
    
test_that("Matrix as input to wmwTest", {
    expect_equivalent(tsMat, tsMatR)
    expect_equivalent(tsMatGreater, tsMatGreaterR)
    expect_equivalent(tsMatLess, tsMatLessR)
    expect_equivalent(tsMat, logicTsMat)
    expect_equivalent(tsMatGreater, logicTsMatGreater)
})

intMat <- matrix(as.integer(rnorm(3000, 5, 3)), nrow=1000)
intMatTs <- wmwTest(intMat, rl, valType="p.two.sided")
intMatGreater <- wmwTest(intMat, rl, valType="p.greater")
intMatLess <- wmwTest(intMat, rl, valType="p.less")
intMatTsR <- apply(intMat, 2, function(x) wmwTestInR(x, rl, valType="p.two.sided"))
intMatGreaterR <- apply(intMat, 2, function(x) wmwTestInR(x, rl, valType="p.greater"))
intMatLessR <- apply(intMat, 2, function(x) wmwTestInR(x, rl, valType="p.less"))

test_that("Integer matrix as input", {
    expect_equivalent(intMatTs, intMatTsR)
    expect_equivalent(intMatGreater, intMatGreaterR)
    expect_equivalent(intMatLess, intMatLessR)
})
    
## transformed val types
tsMatGreaterScore <- wmwTest(rmat, which(rl), valType="abs.log10p.greater")
expTsMatGreaterScore <- abs(log10(tsMatGreater))

tsMatLessScore <- wmwTest(rmat, which(rl), valType="log10p.less")
expTsMatLessScore <- log10(tsMatLess)

tsMatScore <- wmwTest(rmat, which(rl), valType="abs.log10p.two.sided")
expTsMatScore <- abs(log10(tsMat))

matQ <- wmwTest(rmat, which(rl), valType="Q")
expQ <-expTsMatScore
isRevQ <- tsMatLess < tsMatGreater
expQ[isRevQ] <- -expQ[isRevQ]
test_that("Transformed valTypes", {
   expect_equivalent(tsMatGreaterScore,expTsMatGreaterScore)
   expect_equivalent(tsMatLessScore,expTsMatLessScore)
   expect_equivalent(tsMatScore, expTsMatScore)
   expect_equivalent(matQ, expQ)
})


## using ExpressionSet
data(sample.ExpressionSet)
testSet <- sample.ExpressionSet
fData(testSet)$GeneSymbol <- paste("GENE_",1:nrow(testSet), sep="")
mySig1 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.25, 0.75), replace=TRUE)

eSetInt <- wmwTest(testSet, which(mySig1), valType="p.greater")
eSetIntR <- apply(exprs(testSet), 2, function(x)
                  wmwTestInR(x, which(mySig1), valType="p.greater"))



## using lists
mySig2 <- sample(c(TRUE, FALSE), nrow(testSet), prob=c(0.6, 0.4), replace=TRUE)
sigLogLists <- list(first=mySig1, second=mySig2)
sigIntLists <- list(first=which(mySig1), second=which(mySig2))
sigLog <- wmwTest(testSet, sigLogLists, valType="p.greater")
sigInt <- wmwTest(testSet, sigIntLists, valType="p.greater")
sigLogR <- apply(exprs(testSet), 2, function(x)
                 sapply(sigLogLists, function(vec)
                        wmwTestInR(x, vec, valType="p.greater")))

test_that("eSet", {
    expect_equivalent(eSetInt, eSetIntR)
    expect_equivalent(sigLog, sigInt)
    expect_equivalent(sigLog, sigLogR)
})

## using GMT file
gmt_file <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt", package="BioQC")
gmt_list <- readGmt(gmt_file)

gss <- sample(unlist(sapply(gmt_list, function(x) x$genes)), 500)
eset<-new("ExpressionSet",
          exprs=matrix(rnorm(4000), nrow=500L),
          phenoData=new("AnnotatedDataFrame", data.frame(Sample=LETTERS[1:8])),
          featureData=new("AnnotatedDataFrame",data.frame(GeneSymbol=gss)))
esetWmwRes <- wmwTest(eset ,gmt_list, valType="p.greater")
gmtInd <- matchGenes(gmt_list, eset, "GeneSymbol")
esetWmwResR <- wmwTest(exprs(eset), gmtInd, valType="p.greater")
                             
test_that("eSet and GMT", {
    expect_equivalent(esetWmwRes, esetWmwResR)
})

context("Test wmwTest for SignedGenesets")
testGenes <- c("AKT1", "AKT2", "EGFR", "ERBB2", "TSC1", "TSC2", "F4")
testRows <- c(testGenes, paste("Gene", (length(testGenes)+1):100, sep=""))
testRawSignedGenesets <- list(GS1=list(name="GS1",
                                  pos=c("AKT1", "AKT2"),
                                  neg=c("TSC1", "TSC2", "PR")),
                              GS2=list(name="GS2",
                                  pos=c("EGFR","ERBB3"),
                                  neg=c("ERBB2", "ERBB4")),
                              GS3=list(name="GS3",
                                  pos=NULL,
                                  neg=c("TSC1", "TSC2")),
                              GS4=list(name="GS4",
                                  pos=c("EGFR", "ERBB2", "PR", NA),
                                  neg=NULL),
                              GS5=list(name="GS5", pos=NULL, neg=NULL))
testSignedGenesets <- SignedGenesets(testRawSignedGenesets)
testSignedMatch <- matchGenes(testSignedGenesets, testRows)
expSignedMatch <- list(list(pos=1:2, neg=5:6),
                       list(pos=3L, neg=4L),
                       list(pos=NULL, neg=5:6),
                       list(pos=3:4, neg=NULL),
                       list(pos=NULL, neg=NULL))

testMatrix <- matrix(rnorm(1000, sd=0.1),
                     nrow=100,
                     dimnames=list(testRows, NULL))
testMatrix[1,] <- testMatrix[1,]+20
testMatrix[2,] <- testMatrix[2,]+10
testMatrix[3,] <- testMatrix[3,]-40
testMatrix[4,] <- testMatrix[4,]-30
testMatrix[5,] <- testMatrix[5,]-20
testMatrix[6,] <- testMatrix[6,]-10
testMatrixRankHead <- apply(testMatrix, 2, rank)[1:6,]
expMatrixRankHead <- matrix(rep(c(100, 99, 1, 2, 3, 4), each=10), ncol=10, byrow=TRUE)


testSignedGreater <- wmwTest(testMatrix, testSignedMatch, valType="p.greater")
testSignedLess <- wmwTest(testMatrix, testSignedMatch, valType="p.less")
testSignedTwoSided <- wmwTest(testMatrix, testSignedMatch, valType="p.two.sided")
testSignedQ <- wmwTest(testMatrix, testSignedMatch, valType="Q")
testSignedU <- wmwTest(testMatrix, testSignedMatch, valType="U")
## U statistic is calculated by hand
expSignedUBase <- c(0, 99, 4, 196, 0)
expSignedU <- matrix(rep(expSignedUBase, each=10L), ncol=10, byrow=TRUE)

test_that("wmwTest works for signed genesets", {
              expect_equivalent(testMatrixRankHead, expMatrixRankHead)
              expect_equivalent(testSignedMatch, expSignedMatch)
              expect_equivalent(testSignedU, expSignedU)
          })

testSignedUcol1 <- wmwTest(testMatrix[,1], testSignedMatch, valType="U")
testSignedEset <- new("ExpressionSet",
                      exprs=testMatrix)
testSignedEsetU <- wmwTest(testSignedEset, testSignedMatch, valType="U")

test_that("wmwTest works for signed genesets and polymorphism", {
              expect_equivalent(testSignedUcol1, expSignedUBase)
              expect_equivalent(testSignedEsetU, expSignedU)

              expect_equal(names(testSignedUcol1), names(testRawSignedGenesets))
              expect_equal(rownames(testSignedEsetU), names(testRawSignedGenesets))
          })
