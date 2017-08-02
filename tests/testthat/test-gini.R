library(ineq)

vec1 <- c(2.0, 3.0, 4.0, 5.0)
vec2 <- rep(2.0, 4)
vec3 <- c(4.0, 5.0, 2.0, 3.0, NA, 4.4)
vec4 <- c(1L, 2L, 4L, 3L, 5L, NA)
vec5 <- c(0, 1, 2, -1, 3, 4, NA)

context("Test gini for vectors")

test_that("BioQC::gini works as well as ineq::Gini",
          expect_equal(ineq::Gini(vec1),
                       BioQC::gini(vec1)))

test_that("BioQC::gini works as well as ineq::Gini for all-equal cases",
          expect_equal(ineq::Gini(vec2),
                       BioQC::gini(vec2)))

test_that("BioQC::gini works well for unsorted numeric vectors with NA values",
          expect_equal(ineq::Gini(vec3[!is.na(vec3)]),
                       BioQC::gini(vec3)))

test_that("BioQC::gini works well for unsorted integer vectors with NA values",
          expect_equal(ineq::Gini(vec4[!is.na(vec4)]),
                       BioQC::gini(vec4)))

test_that("BioQC::gini reports error if there are negative values",
          expect_error(BioQC::gini(vec5)))

context("Test gini for matrix")

testMat <- rbind(vec1, vec2)
testMatGini <- gini(testMat)
testMatExpGini <- c(ineq::Gini(vec1),
                    ineq::Gini(vec2))
test_that("BioQC::gini works with Matrix",
          expect_equal(testMatGini, testMatExpGini))

context("Test gini for large matrix")
set.seed(1887)
testMat2 <- abs(matrix(rnorm(1000), nrow=10))
testMat2Gini <- gini(testMat2)
testMat2ExpGini <- apply(testMat2, 1, ineq::Gini)
test_that("BioQC::gini works with large matrix",
          expect_equal(testMat2Gini, testMat2ExpGini))
