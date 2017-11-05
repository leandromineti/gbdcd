context("Gaussian BDCD")

context("RcppArma functions")

test_that("Frequency matrix", {
  t <- c(1:64)
  dt <- RcppFreqMatrix(t)
  expect_is(dt, "matrix")
  expect_equal(nrow(dt), 64)
  expect_equal(ncol(dt), 64)
})

data("aneeldata", package = "gbdcd")
data("aneelshape", package = "gbdcd")

test_that("Partitioning", {
  dt <- aneeldata[[1]]
  part <- RcppPartition(dt, c(3, 40)) 
  expect_equal(length(part), 57)
})

context("Gbdcd")

test_that("gbdcd function", {
  res <- gaussianBDCD(y = aneelshape$z_Precipitation, neigh = aneeldata$connections,
                      n_iterations = 1000, burn_in = 500, c = 0.35, 0, sigma0 = sqrt(2))
  expect_is(res, "list")
})
