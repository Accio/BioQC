library(ineq)

vec1 <- c(2.0, 3.0, 4.0, 5.0)
vec2 <- rep(2.0, 4)

context("Test gini")

test_that("BioQC::gini works as well as ineq::Gini",
          expect_equal(ineq::Gini(vec1),
                       BioQC::gini(vec1)))

test_that("BioQC::gini works as well as ineq::Gini for all-equal cases",
          expect_equal(ineq::Gini(vec2),
                       BioQC::gini(vec2)))

context("Test gini for matrix")

testMat <- rbind(vec1, vec2)
testMatGini <- gini(testMat)
testMatExpGini <- c(ineq::Gini(vec1),
                    ineq::Gini(vec2))
test_that("BioQC::gini works with Matrix",
          expect_equal(testMatGini, testMatExpGini))

context("Test gini for large matrix")
set.seed(1887)
testMat2 <- matrix(rnorm(1000), nrow=10)
testMat2Gini <- gini(testMat2)
testMat2ExpGini <- apply(testMat2, 1, ineq::Gini)
test_that("BioQC::gini works with large matrix",
          expect_equal(testMat2Gini, testMat2ExpGini))
