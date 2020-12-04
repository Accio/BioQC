# test rank-biserial correlation and common language effect size
library(BioQC)

context("Test correspondance between U and rank-biserial correlation")
n1 <- 300
n2 <- 200
rd <- c(rnorm(n1), rnorm(n2, -0.2))
ind <- c(rep(FALSE, n1), rep(TRUE, n2))

# one-dimensional vector
wmw_U <- wmwTest(rd, ind, valType = "U")
wmw_U_f <- wmw_U / n1 / n2
wmw_U_r <- 2 * wmw_U / n1 / n2 - 1
wmw_r <- wmwTest(rd, ind, valType = "r")
wmw_f <- wmwTest(rd, ind, valType = "f")


## matrix form
rmat <- matrix(c(rd, rd + rnorm(length(rd)), rd + rnorm(length(rd))), ncol = 3, byrow = FALSE)

wmw_mat_U <- wmwTest(rmat, ind, valType = "U")
wmw_mat_U_f <- wmw_mat_U / n1 / n2
wmw_mat_U_r <- 2 * wmw_mat_U / n1 / n2 - 1
wmw_mat_r <- wmwTest(rmat, ind, valType = "r")
wmw_mat_f <- wmwTest(rmat, ind, valType = "f")

test_that("Common language effect size and rank-biserial correlation values are compatible", {
  expect_equivalent(wmw_r, 2 * wmw_f - 1)
  expect_equivalent(wmw_mat_r, 2 * wmw_mat_f - 1)
})

test_that("Common language effect size corresponds to U value", {
  expect_equivalent(wmw_f, wmw_U_f)
  expect_equivalent(wmw_mat_f, wmw_mat_U_f)
})

test_that("Rank-biserial correlation corresponds to U value", {
  expect_equivalent(wmw_r, wmw_U_r)
  expect_equivalent(wmw_mat_r, wmw_mat_U_r)
})
