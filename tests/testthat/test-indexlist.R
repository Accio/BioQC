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
          })
test_that("test IndexList: offset 0", {
              offset(testIndexList) <- 0L
              expect_equivalent(testIndexList@.Data,
                                list(0:3,
                                     1:4,
                                     NULL,
                                     c(0L, 2L, 4L, 3L),
                                     c(1L, 3L)))
          })
test_that("test IndexList: offset 2", {
              offset(testIndexList) <- 2L
              expect_equivalent(testIndexList@.Data,
                                list(2:5,
                                     3:6,
                                     NULL,
                                     c(2L, 4L, 6L, 5L),
                                     c(3L, 5L)))
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
resSil <- SignedIndexList(inputSil)
test_that("SignedIndexList from list", {
    expect_equal(names(resSil), paste("GS", LETTERS[1:4], sep="_"))
    expect_equal(resSil@keepNA, FALSE)
    expect_equal(resSil@keepDup, FALSE)
    expect_equal(offset(resSil), 1L)
    expect_equivalent(resSil@.Data, list(list(pos=1:3, neg=NULL),
                                         list(pos=1:3, neg=4:7),
                                         list(pos=NULL, neg=c(3L, 5L)),
                                         list(pos=c(1L, 3L, 4L), neg=c(3L, 5L))))
})

test_that("SignedIndexList: offset 0", {
    offset(resSil) <- 0L
    expect_equivalent(resSil@.Data, list(list(pos=0:2, neg=NULL),
                                         list(pos=0:2, neg=3:6),
                                         list(pos=NULL, neg=c(2L, 4L)),
                                         list(pos=c(0L, 2L, 3L), neg=c(2L, 4L))))
})

test_that("SignedIndexList: offset 2", {
    offset(resSil) <- 2L
    expect_equivalent(resSil@.Data, list(list(pos=2:4, neg=NULL),
                                         list(pos=2:4, neg=5:8),
                                         list(pos=NULL, neg=c(4L, 6L)),
                                         list(pos=c(2L, 4L, 5L), neg=c(4L, 6L))))
})

context("Test matchGene for SignedGenesets")
