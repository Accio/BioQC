# test U1 and U2 values
context("Test correspondance between U, U1, U2 and n1 * n2")

set.seed(10)
n1 <- sample(100, 1)
n2 <- sample(100, 1)
x <- rnorm(n1)
y <- rnorm(n2)

W_U <- wilcox.test(x, y)$statistic
W_U_inv <- wilcox.test(y, x)$statistic


U <- wmwTest(c(x, y), seq_along(x), valType = "U")
U1 <- wmwTest(c(x, y), seq_along(x), valType = "U1")
U2 <- wmwTest(c(x, y), seq_along(x), valType = "U2")

U_inv <- wmwTest(c(y, x), seq_along(y), valType = "U")
U1_inv <- wmwTest(c(y, x), seq_along(y), valType = "U1")
U2_inv <- wmwTest(c(y, x), seq_along(y), valType = "U2")


## matrix form
set.seed(11)
ind <- c(rep(FALSE, n1), rep(TRUE, n2))
nc <- 10
rmat <- matrix(rnorm(length(ind) * nc), ncol = nc)

mat_W_U <- apply(rmat, 2, function(v) wilcox.test(v[ind], v[!ind])$statistic)
mat_W_U_inv <- apply(rmat, 2, function(v) wilcox.test(v[!ind], v[ind])$statistic)

mat_U <- wmwTest(rmat, ind, valType = "U")
mat_U1 <- wmwTest(rmat, ind, valType = "U1")
mat_U2 <- wmwTest(rmat, ind, valType = "U2")

mat_U_inv <- wmwTest(rmat, !ind, valType = "U")
mat_U1_inv <- wmwTest(rmat, !ind, valType = "U1")
mat_U2_inv <- wmwTest(rmat, !ind, valType = "U2")

test_that("Test that U equals to U1", {
  expect_equal(U, U1)
  expect_equal(U_inv, U1_inv)
  expect_equal(mat_W_U, mat_U1)
  expect_equal(mat_W_U_inv, mat_U1_inv)
})

test_that("Test that U1 + U1 equals n1 * n2", {
  expect_equal(U1 + U2, n1 * n2)
  expect_equal(U1_inv + U2_inv, n1 * n2)
  expect_equal(mat_U1 + mat_U2, rep(n1 * n2, nc))
  expect_equal(mat_U1_inv + mat_U2_inv, rep(n1 * n2, nc))
})

test_that("Test that U1 corresponds to U2 if x and y are exchanged", {
  expect_equal(U1, U2_inv)
  expect_equal(U2, U1_inv)
  expect_equal(mat_U1, mat_U2_inv)
  expect_equal(mat_U2, mat_U1_inv)
})

test_that("Test that W from the Wilcoxon test equals to U1", {
  expect_equivalent(U1, W_U)
  expect_equivalent(U1_inv, W_U_inv)
  expect_equivalent(mat_U1, mat_W_U)
  expect_equivalent(mat_U1_inv, mat_W_U_inv)
})
