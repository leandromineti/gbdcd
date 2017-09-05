context("Gaussian BDCD")

context("RcppArma functions")

test_that("Get trip summary by Id: all columns", {
  t <- c(1:64)
  dt <- RcppFreqMatrix(t)
  expect_is(dt, "matrix")
  expect_equal(nrow(dt), 64)
  expect_equal(ncol(dt), 64)
})
