## test wmwTest

NROW <- 100
NCOL <- 40
GSCOUNT <- 20
TOL <- 1E-8
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

expect_equal(Cless, Rless, tolerance=TOL)
expect_equal(Cgreater, Rgreater,tolerance=TOL)
expect_equal(object=Ctwosided, expected=Rtwosided, tolerance=TOL, scale=1L)
