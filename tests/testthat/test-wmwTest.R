## test wmwTest

NROW <- 100
NCOL <- 40
GSCOUNT <- 20
GSSIZE <- sapply(1:GSCOUNT, function(x) sample(1:NROW/2, replace=TRUE))
set.seed(1887)
ind <- lapply(1:GSCOUNT, function(i) sample(1:NROW, GSSIZE[i]))
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

expect_equal(Cless, Rless, tolerance=1E-8)
expect_equal(Cgreater, Rgreater,tolerance=1E-8)
## note that sometimes this fails when the p value of R is 1 whiel of C is around 0.996
expect_equal(Ctwosided, Rtwosided, tolerance=1E-8)


## wmwTest(exprs[ind[[19]],16],exprs[,16]) ## reports segmentation error - this shuold not happen (some sanity check is missing)
