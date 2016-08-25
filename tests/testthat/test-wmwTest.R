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
