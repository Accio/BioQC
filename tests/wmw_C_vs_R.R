library(BioQC)
NROW <- 2000
NCOL <- 20
GSCOUNT <- 5
GSSIZE <- sapply(1:GSCOUNT, function(x) sample(1:NROW/2, replace=TRUE))
set.seed(1887)
ind <- lapply(1:GSCOUNT, function(i) sample(1:NROW, GSSIZE[i]))
exprs <- matrix(round(rnorm(NROW*NCOL),4), nrow=NROW)
gc()
system.time(Cres <- wmwTest(exprs, ind, alternative="less"))
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

testWMW <- function(ind.list, matrix, alternative=c("less", "greater", "two.sided")) {
  alternative <- match.arg(alternative)
  system.time(Cres <- wmwTest(matrix, ind.list, alternative=alternative))
  system.time(Rres <- apply(matrix, 2, function(x)
                            sapply(ind.list, function(y)
                                   wmwTestR(x, y, alternative=alternative, stat=FALSE))))
  if(!all(abs(Cres - Rres)<1E-5)) {
    cat("ERROR\n")
    print(Cres)
    print(Rres)
    print(summary(as.vector(Cres-Rres)))
    quit(save="no", status=1)
  }
}

testWMW(ind, exprs, "less")
testWMW(ind, exprs, "greater")
testWMW(ind, exprs, "two.sided")
