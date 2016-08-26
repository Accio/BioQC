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

testLogical <- c(FALSE, TRUE, FALSE, TRUE)
testLogIndexList <- IndexList(testLogical)
testIntIndexList <- IndexList(which(testLogical))
test_that("test IndexList from logical vector", {
    expect_equivalent(testLogIndexList@.Data, list(c(2,4)))
    expect_equivalent(testIntIndexList@.Data, list(c(2,4)))
})

context("Test SignedIndexList")

inputSil <- list("GS_A"=list(pos=1:3, neg=NULL),
                 "GS_B"=list(pos=1:3, neg=4:7),
                 "GS_C"=list(pos=NULL, neg=c(3L,5L,NA)),
                 "GS_D"=list(pos=c(1L,3L,3L,4L), neg=c(3L,5L,NA)))
## resSil <- SignedIndexList(inputSil)
