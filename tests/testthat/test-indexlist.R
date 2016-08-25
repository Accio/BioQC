library(BioQC)

context("Test IndexList")

testList <- list(GS_A=c(1,2,3,4),
                 GS_B=c(2,3,4,5),
                 GS_C=NULL,
                 GS_D=c(1,3,5,NA),
                 GS_E=c(2,4))
testIndexList <- IndexList(testList, offset=1L)

test_that("test IndexList", {
              expect_equivalent(testList, testIndexList@.Data)
              expect_equal(offset(testIndexList), 1L)
              offset(testIndexList) <- 0L
              expect_equivalent(sapply(testList, function(x) x-1L), testIndexList@.Data)
              offset(testIndexList) <- 2L
              expect_equivalent(sapply(testList, function(x) x+1L), testIndexList@.Data)
          })
