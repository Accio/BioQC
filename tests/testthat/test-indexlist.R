library(BioQC)

context("Test IndexList")

testList <- list(GS_A=c(1,2,3,4,3),
                 GS_B=c(2,3,4,5),
                 GS_C=NULL,
                 GS_D=c(1,3,5,NA,4,5),
                 GS_E=c(2,4))
expList <- list(GS_A=c(1,2,3,4),
                GS_B=c(2,3,4,5),
                GS_C=NULL,
                GS_D=c(1,3,5,4),
                GS_E=c(2,4))
expListNA <- list(GS_A=c(1,2,3,4),
                  GS_B=c(2,3,4,5),
                  GS_C=NULL,
                  GS_D=c(1,3,5, NA ,4),
                  GS_E=c(2,4))
testIndexList <- IndexList(testList, offset=1L)
testIndexListNA <- IndexList(testList, offset=1L, keepNA=TRUE)

test_that("test IndexList", {
              expect_equivalent(expList, testIndexList@.Data)
              expect_equivalent(expListNA, testIndexListNA@.Data)
              expect_equal(offset(testIndexList), 1L)
              offset(testIndexList) <- 0L
              expect_equivalent(sapply(expList, function(x) x-1L), testIndexList@.Data)
              offset(testIndexList) <- 2L
              expect_equivalent(sapply(expList, function(x) x+1L), testIndexList@.Data)
          })
