library(BioQC)
library(testthat)

context("Testing subetting of GmtList")

myList <- GmtList(list("A"=LETTERS[1:3],
            "B"=letters[1:3],
            "C"=LETTERS[5:8]))

test_that("Subsetting of a GmtList with [", {
  expect_identical(myList["A"],
                   myList[1])
  expect_identical(myList[c("A", "C")],
                   myList[c(1,3)])
  expect_warning(myList[c("A", "C", "D")])
})

test_that("Subsetting of a GmtList with [[", {
  expect_identical(myList[["A"]],
                   list(name="A", desc=NULL,
                        genes=LETTERS[1:3], namespace=NULL))
  expect_warning(expect_error(myList[["D"]]))
})