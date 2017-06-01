library(DEXICA)
library(fastICA)
context("Test module prediction")

# Prepare test data
set.seed(1)
x = matrix(rnorm(100), 10, 10)
x = preprocessMatrix(x)
set.seed(1)
m = predictModules(x, n.comp = 3, w.init = 1)
set.seed(1)
n = predictModules(x, n.comp = 3, w.init = 1)
set.seed(1)
o = fastICA(x, n.comp = 3, verbose = FALSE)
set.seed(1)
p = predictModules(x, n.comp = 3, alg.typ = "def")
set.seed(1)
q = fastICA(x, n.comp = 3, alg.typ = "deflation", verbose = FALSE)

test_that("the dimensions of S and A are as expected", {
  expect_equal(ncol(m$S), 3)
  expect_equal(ncol(m$A), 10)
  expect_equal(nrow(m$S), 10)
  expect_equal(nrow(m$A), 3)
})

test_that("S has been normalized", {
  expect_equal(ncol(m$S), 3)
  expect_equal(ncol(m$A), 10)
  expect_equal(nrow(m$S), 10)
  expect_equal(nrow(m$A), 3)
})

test_that("fastICA produces the same result", {
  expect_true(all(round(abs(t(n$S) %*% o$S)) %in% c(10,0)))
  expect_true(all(round(abs(t(p$S) %*% q$S)) %in% c(10,0)))
})

test_that("two runs with the same seed produce the same result", {
  expect_identical(m, n)
})
