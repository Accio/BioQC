## test wmwTest
library(BioQC)

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
system.time(Cless <- wmwTest(exprs, ind, alternative="less"))
system.time(Cgreater <- wmwTest(exprs, ind, alternative="greater"))
system.time(Ctwosided <- wmwTest(exprs, ind, alternative="two.sided"))
wmwTestR <- function(matrix, index, alternative, stat) {
  sub <- rep(FALSE, length(matrix))
  sub[index]=TRUE
  BioQC:::wmw.test(matrix, sub, alternative, stat)
}
system.time(Rless <- apply(exprs, 2, function(x)
                          sapply(ind, function(y) 
                                 wmwTestR(x, y, alternative="less", stat=FALSE))))
system.time(Rgreater <- apply(exprs, 2, function(x)
                              sapply(ind, function(y)
                                     wmwTestR(x, y, alternative="greater", stat=FALSE))))
system.time(Rtwosided <- apply(exprs, 2, function(x)
                              sapply(ind, function(y)
                                     wmwTestR(x, y, alternative="two.sided", stat=FALSE))))

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
