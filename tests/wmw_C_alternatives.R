library(BioQC)
NROW <- 2000
NCOL <- 20
GSCOUNT <- 5
GSSIZE <- sapply(1:GSCOUNT, function(x) sample(1:NROW/2, replace=TRUE))
set.seed(1887)
ind <- lapply(1:GSCOUNT, function(i) sample(1:NROW, GSSIZE[i]))
exprs <- matrix(round(rnorm(NROW*NCOL),4), nrow=NROW)

assertEqual <- function(x,y) stopifnot(identical(x,y))
system.time(Cless <- wmwTest(exprs, ind, alternative="less"))
system.time(ClessL <- wmwTest(exprs, ind, alternative="log10.less"))
assertEqual(log10(Cless), ClessL)

system.time(Cgreater <- wmwTest(exprs, ind, alternative="greater"))
system.time(CgreaterL <- wmwTest(exprs, ind, alternative="abs.log10.greater"))
assertEqual(abs(log10(Cgreater)), CgreaterL)

system.time(Cts <- wmwTest(exprs, ind, alternative="two.sided"))
system.time(CtsL <- wmwTest(exprs, ind, alternative="abs.log10.two.sided"))
assertEqual(abs(log10(Cts)), CtsL)

system.time(Cu <- wmwTest(exprs, ind, alternative="U"))
system.time(Cq <- wmwTest(exprs, ind, alternative="Q"))

manualQ <- matrix(NA, nrow=nrow(Cgreater), ncol=ncol(Cgreater))
isLess <- Cless < Cgreater
manualQ[isLess] <- log10(Cts[isLess])
manualQ[!isLess] <- abs(log10(Cts[!isLess]))
assertEqual(Cq, manualQ)
