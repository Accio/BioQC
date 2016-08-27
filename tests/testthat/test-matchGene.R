context("Test matchGene")

testGenes <- c("AKT1", "AKT2", "EGFR", "ERBB2", "TSC1", "TSC2", NA)
testRawGmtList <- list(GS1=c("AKT1", "AKT2"),
                       GS2=c("EGFR", "ERBB2", "ERBB3", "ERBB4"),
                       GS3=NULL,
                       GS4=c("TSC1", "TSC2", "PR", NA))
testGmtList <- as.gmtlist(testRawGmtList)

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
