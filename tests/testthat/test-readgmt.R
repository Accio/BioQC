library(BioQC)

## readGmt
context("Read in gmt file into a gmtlist object")
testFile <- system.file("extdata/test.gmt", package="BioQC")
testGmt <- readGmt(testFile)
expGmt <- list(list(name="GS_A", desc="TestGeneSet_A", genes=c("AKT1", "AKT2", "AKT3")),
               list(name="GS_B", desc="TestGeneSet_B", genes=c("MAPK1", "MAPK3", "MAPK8")),
               list(name="GS_C_UP", desc="TestGeneSet_C", genes=c("ERBB2", "ERBB3")),
               list(name="GS_C_DN", desc="TestGeneSet_C", genes=c("EGFR", "ERBB4")),
               list(name="GS_D_UP", desc="TestGeneSet_D", genes=c("GATA2", "GATA4")),
               list(name="GS_E_DN", desc="TestGeneSet_E", genes=c("TSC1", "TSC2")))
test_that("readGmt",{
             expect_equivalent(expGmt, testGmt)
         })

## gmtlist2signedGenesets
context("Convert gmtlist (a memory-copy of a GMT file) to a list of signed gene sets")

testInputList <- list(list(name="GeneSetA_UP",genes=LETTERS[1:3]),
                      list(name="GeneSetA_DN", genes=LETTERS[4:6]),
                      list(name="GeneSetB", genes=LETTERS[2:4]),
                      list(name="GeneSetC_DN", genes=LETTERS[1:3]),
                      list(name="GeneSetD_UP", genes=LETTERS[1:3]))
outList.ignore <- gmtlist2signedGenesets(testInputList, nomatch="ignore")
outList.pos <- gmtlist2signedGenesets(testInputList, nomatch="pos")
outList.neg <- gmtlist2signedGenesets(testInputList, nomatch="neg")

exp.ignore <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                   "GeneSetB"=list(pos=NULL, neg=NULL),
                   "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                   "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

exp.pos <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                "GeneSetB"=list(pos=LETTERS[2:4], neg=NULL),
                "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

exp.neg <- list("GeneSetA"=list(pos=LETTERS[1:3], neg=LETTERS[4:6]),
                "GeneSetB"=list(pos=NULL, neg=LETTERS[2:4]),
                "GeneSetC"=list(pos=NULL, neg=LETTERS[1:3]),
                "GeneSetD"=list(pos=LETTERS[1:3], neg=NULL))

test_that("gmtlist2signedGenesets, ignore non-matching genesets", {
              expect_equivalent(outList.ignore, exp.ignore)
          })

test_that("gmtlist2signedGenesets, ignore non-matching as positive", {
             expect_equivalent(outList.pos, exp.pos)
         })
test_that("gmtlist2signedGenesets, ignore non-matching as negative", {
             expect_equivalent(outList.neg, exp.neg)
         })

## readSignedGmt
context("Read in gmt file into a signed_genesets object")
testSignedGenesets.ignore <- readSignedGmt(testFile, nomatch="ignore")
testSignedGenesets.pos <- readSignedGmt(testFile, nomatch="pos")
testSignedGenesets.neg <- readSignedGmt(testFile, nomatch="neg")

expSignedGenesets.ignore <- list(GS_A=list(pos=NULL, neg=NULL),
                                 GS_B=list(pos=NULL, neg=NULL),
                                 GS_C=list(pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4")),
                                 GS_D=list(pos=c("GATA2", "GATA4"), neg=NULL),
                                 GS_E=list(pos=NULL, neg=c("TSC1", "TSC2")))
expSignedGenesets.pos <- list(GS_A=list(pos=c("AKT1", "AKT2", "AKT3"), neg=NULL),
                              GS_B=list(pos=c("MAPK1", "MAPK3", "MAPK8"), neg=NULL),
                              GS_C=list(pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4")),
                              GS_D=list(pos=c("GATA2", "GATA4"), neg=NULL),
                              GS_E=list(pos=NULL, neg=c("TSC1", "TSC2")))
expSignedGenesets.neg <- list(GS_A=list(pos=NULL, neg=c("AKT1", "AKT2", "AKT3")),
                              GS_B=list(pos=NULL, neg=c("MAPK1", "MAPK3", "MAPK8")),
                              GS_C=list(pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4")),
                              GS_D=list(pos=c("GATA2", "GATA4"), neg=NULL),
                              GS_E=list(pos=NULL, neg=c("TSC1", "TSC2")))
test_that("readSignedGmt, nomatch ingore",{
              expect_equivalent(testSignedGenesets.ignore, expSignedGenesets.ignore)
          })
test_that("readSignedGmt, nomatch pos",{
              expect_equivalent(testSignedGenesets.pos, expSignedGenesets.pos)
          })
test_that("readSignedGmt, nomatch neg",{
              expect_equivalent(testSignedGenesets.neg, expSignedGenesets.neg)
          })
