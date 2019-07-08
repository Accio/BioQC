library(BioQC)

## readGmt
context("Read gmt file into a GmtList object")
testFile <- system.file("extdata/test.gmt", package="BioQC")
testGmt <- readGmt(testFile)
expGmt <- list(list(name="GS_A", desc="TestGeneSet_A", genes=c("AKT1", "AKT2", "AKT3"), namespace=NULL),
               list(name="GS_B", desc="TestGeneSet_B", genes=c("MAPK1", "MAPK3", "MAPK8"), namespace=NULL),
               list(name="GS_C_UP", desc="TestGeneSet_C", genes=c("ERBB2", "ERBB3"), namespace=NULL),
               list(name="GS_C_DN", desc="TestGeneSet_C", genes=c("EGFR", "ERBB4"), namespace=NULL),
               list(name="GS_D_UP", desc="TestGeneSet_D", genes=c("GATA2", "GATA4"), namespace=NULL),
               list(name="GS_D_DN", desc="TestGeneSet_D", genes=c("GATA1", "GATA3"), namespace=NULL),
               list(name="GS_E_DN", desc="TestGeneSet_E", genes=c("TSC1", "TSC2"), namespace=NULL))
test_that("readGmt reads one GMT file without namespace",{
             expect_equivalent(expGmt, testGmt)
         })

context("Read gmt file into a GmtList object with namespace")
testFile <- system.file("extdata/test.gmt", package="BioQC")
testGmt <- readGmt(testFile, namespace="BioQC")
expGmt <- list(list(name="GS_A", desc="TestGeneSet_A", genes=c("AKT1", "AKT2", "AKT3"), namespace="BioQC"),
               list(name="GS_B", desc="TestGeneSet_B", genes=c("MAPK1", "MAPK3", "MAPK8"), namespace="BioQC"),
               list(name="GS_C_UP", desc="TestGeneSet_C", genes=c("ERBB2", "ERBB3"), namespace="BioQC"),
               list(name="GS_C_DN", desc="TestGeneSet_C", genes=c("EGFR", "ERBB4"), namespace="BioQC"),
               list(name="GS_D_UP", desc="TestGeneSet_D", genes=c("GATA2", "GATA4"), namespace="BioQC"),
               list(name="GS_D_DN", desc="TestGeneSet_D", genes=c("GATA1", "GATA3"), namespace="BioQC"),
               list(name="GS_E_DN", desc="TestGeneSet_E", genes=c("TSC1", "TSC2"), namespace="BioQC"))
test_that("readGmt reads one GMT file with namespace",{
  expect_equivalent(expGmt, testGmt)
})

context("Read two gmt files into a GmtList object with namespaces")
testFile <- system.file("extdata/test.gmt", package="BioQC")
testGmtNamed <- readGmt(ns1=testFile, ns2=testFile)
testGmtUnnamed <- readGmt(testFile, testFile, namespace=c("ns1", "ns2"))
expGmtTwoNs <- list(list(name="GS_A", desc="TestGeneSet_A", genes=c("AKT1", "AKT2", "AKT3"), namespace="ns1"),
               list(name="GS_B", desc="TestGeneSet_B", genes=c("MAPK1", "MAPK3", "MAPK8"), namespace="ns1"),
               list(name="GS_C_UP", desc="TestGeneSet_C", genes=c("ERBB2", "ERBB3"), namespace="ns1"),
               list(name="GS_C_DN", desc="TestGeneSet_C", genes=c("EGFR", "ERBB4"), namespace="ns1"),
               list(name="GS_D_UP", desc="TestGeneSet_D", genes=c("GATA2", "GATA4"), namespace="ns1"),
               list(name="GS_D_DN", desc="TestGeneSet_D", genes=c("GATA1", "GATA3"), namespace="ns1"),
               list(name="GS_E_DN", desc="TestGeneSet_E", genes=c("TSC1", "TSC2"), namespace="ns1"),
               list(name="GS_A", desc="TestGeneSet_A", genes=c("AKT1", "AKT2", "AKT3"), namespace="ns2"),
               list(name="GS_B", desc="TestGeneSet_B", genes=c("MAPK1", "MAPK3", "MAPK8"), namespace="ns2"),
               list(name="GS_C_UP", desc="TestGeneSet_C", genes=c("ERBB2", "ERBB3"), namespace="ns2"),
               list(name="GS_C_DN", desc="TestGeneSet_C", genes=c("EGFR", "ERBB4"), namespace="ns2"),
               list(name="GS_D_UP", desc="TestGeneSet_D", genes=c("GATA2", "GATA4"), namespace="ns2"),
               list(name="GS_D_DN", desc="TestGeneSet_D", genes=c("GATA1", "GATA3"), namespace="ns2"),
               list(name="GS_E_DN", desc="TestGeneSet_E", genes=c("TSC1", "TSC2"), namespace="ns2"))
test_that("readGmt reads two GMT files with namespace",{
  expect_equivalent(expGmtTwoNs, testGmtNamed)
  expect_equivalent(expGmtTwoNs, testGmtUnnamed)
})


## appendGmtList
context("appendGmtList works as expected")
testFile <- system.file("extdata/test.gmt", package="BioQC")
testGmt1 <- readGmt(testFile, namespace="ns1")
testGmt2 <- readGmt(testFile, namespace = "ns2")
testGmt3 <- readGmt(testFile, namespace = "ns3")
testGmtAppended <- appendGmtList(testGmt1, testGmt2, testGmt3)
testGmtAppendedExp <- readGmt(ns1=testFile, ns2=testFile, ns3=testFile)
test_that("appendGmtList works as expected for three GmtList objects",{
  expect_equivalent(testGmtAppended, testGmtAppendedExp)
})

context("Read gmt file into a SignedGenesets object")
testSignedGenesets <- readSignedGmt(testFile, nomatch="pos")
expSignedGenesets <- list(list(name="GS_A", pos=c("AKT1", "AKT2", "AKT3"), neg=NULL, namespace=NULL),
                          list(name="GS_B", pos=c("MAPK1", "MAPK3", "MAPK8"), neg=NULL, namespace=NULL),
                          list(name="GS_C", pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4"), namespace=NULL),
                          list(name="GS_D", pos=c("GATA2", "GATA4"), neg=c("GATA1", "GATA3"), namespace=NULL),
                          list(name="GS_E", pos=NULL, neg=c("TSC1", "TSC2"), namespace=NULL))
test_that("readSignedGmt", {
              expect_equal(expSignedGenesets, testSignedGenesets@.Data)
          })

## as.gmtlist
context("Convert a list of gene symbols into a gmtlist object")
testVec <- list(GeneSet1=c("AKT1", "AKT2"),
                GeneSet2=c("MAPK1", "MAPK3"),
                GeneSet3=NULL)
testVecGmtlist <- as.gmtlist(testVec)
expVecGmtlist <- list(GeneSet1=list(name="GeneSet1", desc=NULL, genes=c("AKT1", "AKT2"), namespace=NULL),
                      GeneSet2=list(name="GeneSet2", desc=NULL, genes=c("MAPK1", "MAPK3"), namespace=NULL),
                      GeneSet3=list(name="GeneSet3", desc=NULL, genes=NULL, namespace=NULL))

testVecGmtlist.desc <- as.gmtlist(testVec, desc=c("GS1", "GS2", "GS3"))
expVecGmtlist.desc <- list(GeneSet1=list(name="GeneSet1", desc="GS1", genes=c("AKT1", "AKT2"), namespace=NULL),
                      GeneSet2=list(name="GeneSet2", desc="GS2", genes=c("MAPK1", "MAPK3"), namespace=NULL),
                      GeneSet3=list(name="GeneSet3", desc="GS3", genes=NULL, namespace=NULL))
test_that("as.gmtlist",{
              expect_equivalent(testVecGmtlist, expVecGmtlist)
              expect_equivalent(testVecGmtlist.desc, expVecGmtlist.desc)
         })

## gmtlist2signedGenesets
context("Convert gmtlist (a memory-copy of a GMT file) to a list of signed gene sets")

testInputList <- list(list(name="GeneSetA_UP",genes=LETTERS[1:3], namespace=NULL),
                      list(name="GeneSetA_DN", genes=LETTERS[4:6], namespace=NULL),
                      list(name="GeneSetB", genes=LETTERS[2:4], namespace=NULL),
                      list(name="GeneSetC_DN", genes=LETTERS[1:3], namespace=NULL),
                      list(name="GeneSetD_UP", genes=LETTERS[1:3], namespace=NULL))
outList.ignore <- gmtlist2signedGenesets(testInputList, nomatch="ignore")
outList.pos <- gmtlist2signedGenesets(testInputList, nomatch="pos")
outList.neg <- gmtlist2signedGenesets(testInputList, nomatch="neg")

exp.ignore <- list("GeneSetA"=list(name="GeneSetA", pos=LETTERS[1:3], neg=LETTERS[4:6], namespace=NULL),
                   "GeneSetB"=list(name="GeneSetB", pos=NULL, neg=NULL, namespace=NULL),
                   "GeneSetC"=list(name="GeneSetC", pos=NULL, neg=LETTERS[1:3], namespace=NULL),
                   "GeneSetD"=list(name="GeneSetD", pos=LETTERS[1:3], neg=NULL, namespace=NULL))

exp.pos <- list("GeneSetA"=list(name="GeneSetA", pos=LETTERS[1:3], neg=LETTERS[4:6], namespace=NULL),
                "GeneSetB"=list(name="GeneSetB", pos=LETTERS[2:4], neg=NULL, namespace=NULL),
                "GeneSetC"=list(name="GeneSetC", pos=NULL, neg=LETTERS[1:3], namespace=NULL),
                "GeneSetD"=list(name="GeneSetD", pos=LETTERS[1:3], neg=NULL, namespace=NULL))

exp.neg <- list("GeneSetA"=list(name="GeneSetA", pos=LETTERS[1:3], neg=LETTERS[4:6], namespace=NULL),
                "GeneSetB"=list(name="GeneSetB",pos=NULL, neg=LETTERS[2:4], namespace=NULL),
                "GeneSetC"=list(name="GeneSetC", pos=NULL, neg=LETTERS[1:3], namespace=NULL),
                "GeneSetD"=list(name="GeneSetD", pos=LETTERS[1:3], neg=NULL, namespace=NULL))

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

expSignedGenesets.ignore <- list(GS_A=list(name="GS_A", pos=NULL, neg=NULL, namespace=NULL),
                                 GS_B=list(name="GS_B", pos=NULL, neg=NULL, namespace=NULL),
                                 GS_C=list(name="GS_C", pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4"), namespace=NULL),
                                 GS_D=list(name="GS_D", pos=c("GATA2", "GATA4"), neg=c("GATA1", "GATA3"), namespace=NULL),
                                 GS_E=list(name="GS_E", pos=NULL, neg=c("TSC1", "TSC2"), namespace=NULL))
expSignedGenesets.pos <- list(GS_A=list(name="GS_A", pos=c("AKT1", "AKT2", "AKT3"), neg=NULL, namespace=NULL),
                              GS_B=list(name="GS_B", pos=c("MAPK1", "MAPK3", "MAPK8"), neg=NULL, namespace=NULL),
                              GS_C=list(name="GS_C", pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4"), namespace=NULL),
                              GS_D=list(name="GS_D", pos=c("GATA2", "GATA4"), neg=c("GATA1", "GATA3"), namespace=NULL),
                              GS_E=list(name="GS_E", pos=NULL, neg=c("TSC1", "TSC2"), namespace=NULL))
expSignedGenesets.neg <- list(GS_A=list(name="GS_A", pos=NULL, neg=c("AKT1", "AKT2", "AKT3"), namespace=NULL),
                              GS_B=list(name="GS_B", pos=NULL, neg=c("MAPK1", "MAPK3", "MAPK8"), namespace=NULL),
                              GS_C=list(name="GS_C", pos=c("ERBB2", "ERBB3"), neg=c("EGFR", "ERBB4"), namespace=NULL),
                              GS_D=list(name="GS_D", pos=c("GATA2", "GATA4"), neg=c("GATA1", "GATA3"), namespace=NULL),
                              GS_E=list(name="GS_E", pos=NULL, neg=c("TSC1", "TSC2"), namespace=NULL))
test_that("readSignedGmt, nomatch ingore",{
              expect_equivalent(testSignedGenesets.ignore, expSignedGenesets.ignore)
          })
test_that("readSignedGmt, nomatch pos",{
              expect_equivalent(testSignedGenesets.pos, expSignedGenesets.pos)
          })
test_that("readSignedGmt, nomatch neg",{
              expect_equivalent(testSignedGenesets.neg, expSignedGenesets.neg)
          })
