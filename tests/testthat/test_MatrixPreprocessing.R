library(DEXICA)
context("Test matrix preprocessing")

# Prepare test data
set.seed(1)
x = matrix(rnorm(100), 10, 10)
x.p = preprocessMatrix(x)

test_that("the means of columns of x are 0", {
  expect_true(all(apply(x.p, 2, FUN = function(y) round(mean(y), 6)) == 0))
})

x.p = preprocessMatrix(x, center.cols = TRUE, scale.cols = TRUE)
test_that("the SDs of columns of x are 1", {
  expect_true(all(apply(x.p, 2, FUN = function(y) round(sd(y), 6)) == 1))
})

x.p = preprocessMatrix(x, center.cols = TRUE, scale.cols = TRUE, center.rows = TRUE)
test_that("the means of rows of x are 0", {
  expect_true(all(apply(x.p, 1, FUN = function(y) round(mean(y), 6)) == 0))
})

x.p = preprocessMatrix(x, center.cols = TRUE, scale.cols = TRUE, center.rows = TRUE, scale.rows = TRUE)
test_that("the SDs of rows of x are 1", {
  expect_true(all(apply(x.p, 1, FUN = function(y) round(sd(y), 6)) == 1))
})

x.p = preprocessMatrix(x, center.cols = TRUE, scale.cols = TRUE, center.rows = TRUE,
                       scale.rows = TRUE, simultaneous = TRUE)
test_that("simultaneous row/col scaling works", {
  expect_true(all(apply(x.p, 1, FUN = function(y) round(sd(y), 6)) == 1))
  expect_true(all(apply(x.p, 1, FUN = function(y) round(mean(y), 6)) == 0))
  expect_true(all(apply(x.p, 2, FUN = function(y) round(sd(y), 6)) == 1))
  expect_true(all(apply(x.p, 2, FUN = function(y) round(mean(y), 6)) == 0))
})

