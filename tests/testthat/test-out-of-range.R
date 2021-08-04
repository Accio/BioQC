## test out-of-range errors
mat <- matrix(rnorm(100), nrow=10)

## well-behaving indices
wbids <- list(gs1=c(1,3,5),
	 gs2=c(1),
	 gs3=1:10)
wmwTest(mat, wbids)

test_that("index out of range errors are found", {
  expect_error(wmwTest(mat, -1),
               "Index out of range: gene set 1, gene 1")
  expect_error(wmwTest(mat, 11),
               "Index out of range: gene set 1, gene 1")
  expect_error(wmwTest(mat, list(1, -1)),
               "Index out of range: gene set 2, gene 1")
	       })

