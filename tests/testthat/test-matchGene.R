context("Test matchGene for GmtList")

testGenes <- c("AKT1", "AKT2", "EGFR", "ERBB2", "TSC1", "TSC2", NA)
testCharQuery <- c("AKT1", "AKT2", "AKT3")
expCharQuery <- list(TempGeneSet=c(1L,2L))
test_that("matchGenes works for characters and characters", {
              expect_equivalent(matchGenes(testCharQuery, testGenes), expCharQuery)
          })

testRawGmtList <- list(GS1=c("AKT1", "AKT2"),
                       GS2=c("EGFR", "ERBB2", "ERBB3", "ERBB4"),
                       GS3=NULL,
                       GS4=c("TSC1", "TSC2", "PR", NA))
testGmtList <- as.GmtList(testRawGmtList)
testGmtListWithMethod <- GmtList(testRawGmtList)
test_that("GmtList works as a wrapper of as.GmtList", {
              expect_equal(testGmtList, testGmtListWithMethod)
          })

testMatch <- matchGenes(testGmtList, testGenes)
expMatchList <- list(1:2, 3:4, NULL, 5:6)

test_that("matchGenes works for GmtList and character", {
    expect_equal(offset(testMatch), 1L)
    expect_equivalent(testMatch@.Data, expMatchList)
})

testMatrix <- matrix(rnorm(49),nrow=length(testGenes), dimnames=list(testGenes, NULL))
testMatrixMatch <- matchGenes(testGmtList, testMatrix)
test_that("matchGenes works for GmtList and matrix", {
    expect_equal(offset(testMatrixMatch), 1L)
    expect_equivalent(testMatrixMatch@.Data, expMatchList)
})

testMatrixCharMatch <- matchGenes(testCharQuery, testMatrix)
test_that("matchGenes works for characters and matrix", {
    expect_equal(offset(testMatrixCharMatch), 1L)
    expect_equivalent(testMatrixCharMatch@.Data, expCharQuery)
})


testMatrix2 <- testMatrix; rownames(testMatrix2)[is.na(rownames(testMatrix2))] <- ""
testEset <- new("ExpressionSet",
                exprs=testMatrix2,
                featureData=new("AnnotatedDataFrame",
                    data.frame(GeneSymbol=testGenes,
                               myID=testGenes,
                               row.names=rownames(testMatrix2))))
testEsetMatch.GeneSymbol <- matchGenes(testGmtList, testEset)
testEsetMatch.myID <- matchGenes(testGmtList, testEset, "myID")
testEsetMatch.fname <- matchGenes(testGmtList, testEset, NULL)
test_that("matchGenes works for GmtList and eSet", {
    expect_equal(offset(testEsetMatch.GeneSymbol), 1L)
    expect_equal(offset(testEsetMatch.myID), 1L)
    expect_equal(offset(testEsetMatch.fname), 1L)

    expect_equivalent(testEsetMatch.GeneSymbol@.Data, expMatchList)
    expect_equivalent(testEsetMatch.myID@.Data, expMatchList)
    expect_equivalent(testEsetMatch.fname@.Data, expMatchList)
})

testEsetCharMatch <- matchGenes(testCharQuery, testEset)
test_that("matchGenes works for character and eSet", {
    expect_equivalent(testEsetCharMatch, expCharQuery)
})
context("Test matchGene for SignedGenesets")
testRawSignedGenesets <- list(GS1=list(name="GS1",
                                  pos=c("AKT1", "AKT2"),
                                  neg=c("TSC1", "TSC2", "PR")),
                              GS2=list(name="GS2",
                                  pos=c("EGFR","ERBB3"),
                                  neg=c("ERBB2", "ERBB4")),
                              GS3=list(name="GS3",
                                  pos=NULL,
                                  neg=c("TSC1", "TSC2")),
                              GS4=list(name="GS4",
                                  pos=c("ERBB2", "ERBB4", "PR", NA),
                                  neg=NULL),
                              GS5=list(name="GS5", pos=NULL, neg=NULL))
testSignedGenesets <- SignedGenesets(testRawSignedGenesets)
testSignedMatch <- matchGenes(testSignedGenesets, testGenes)
expSignedMatchList <- list(GS1=list(pos=1:2, neg=5:6),
                           GS2=list(pos=3L, neg=4L),
                           GS3=list(pos=NULL, neg=5:6),
                           GS4=list(pos=4L, neg=NULL),
                           GS5=list(pos=NULL, neg=NULL))

test_that("matchGenes works for SignedGenesets and character", {
    expect_equal(offset(testSignedMatch), 1L)
    expect_equivalent(testSignedMatch@.Data, expSignedMatchList)
    expect_equal(names(testSignedMatch), names(expSignedMatchList))
})

testMatrixSignedMatch <- matchGenes(testSignedGenesets, testMatrix)
test_that("matchSignedGenes works for SignedGenesets and matrix", {
    expect_equal(offset(testMatrixSignedMatch), 1L)
    expect_equivalent(testMatrixSignedMatch@.Data, expSignedMatchList)
})

testEsetSignedMatch.GeneSymbol <- matchGenes(testSignedGenesets, testEset)
testEsetSignedMatch.myID <- matchGenes(testSignedGenesets, testEset, col="myID")
testEsetSignedMatch.fname <- matchGenes(testSignedGenesets, testEset, col=NULL)
test_that("matchSignedGenes works for SignedGenesets and eSet", {
    expect_equal(offset(testEsetSignedMatch.GeneSymbol), 1L)
    expect_equal(offset(testEsetSignedMatch.myID), 1L)
    expect_equal(offset(testEsetSignedMatch.fname), 1L)

    expect_equal(names(testEsetSignedMatch.GeneSymbol), names(expSignedMatchList))
    expect_equal(names(testEsetSignedMatch.myID), names(expSignedMatchList))
    expect_equal(names(testEsetSignedMatch.fname), names(expSignedMatchList))
    
    expect_equivalent(testEsetSignedMatch.GeneSymbol@.Data, expSignedMatchList)
    expect_equivalent(testEsetSignedMatch.myID@.Data, expSignedMatchList)
    expect_equivalent(testEsetSignedMatch.fname@.Data, expSignedMatchList)
})
