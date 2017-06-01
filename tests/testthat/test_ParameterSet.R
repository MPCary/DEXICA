library(DEXICA)
context("Test ParameterSet creation and functionality")

# Simplest case
p = parameterSet()

test_that("the ParamterSet object was created successfully", {
  expect_equal(class(p)[[1]], "ParameterSet")
  expect_equal(getParamSetCount(p), 1)
  expect_error(getOneParamSet(p, 0))
  expect_error(getOneParamSet(p, 2))
  expect_error(getOneParamSet(p, "the_first_one"))
})


# Add some complexity
p = parameterSet(n.comp = seq(from = 5, to = 100, by = 5), center.cols = c(TRUE, FALSE))

test_that("the ParamterSet object has a correct param.space table", {
  expect_equal(getParamSetCount(p), 40)
  expect_equal(length(unique(getAllParamSets(p)[, 'n.comp'])), 20)
  expect_equal(length(unique(getAllParamSets(p)[, 'center.cols'])), 2)
})


# Test w.init
p = parameterSet(w.init = c(-1,0,1))

test_that("w.init is numeric and the proper values", {
  expect_equal(class(getOneParamSet(p, 1)['w.init'][[1]]), "numeric")
  expect_equal(getAllParamSets(p)[, 'w.init'], c(-1,0,1))
})

# Test reseeding
p = reseedParamSets(p)
test_that("w.init is numeric and the proper values", {
  expect_equal(class(getOneParamSet(p, 1)['w.init'][[1]]), "numeric")
  expect_false(all(getAllParamSets(p)[, 'w.init'] %in% c(-1,0,1)))
})

p = reseedParamSets(p, new.seeds = c(4,5,6))
test_that("w.init is numeric and the proper values", {
  expect_equal(class(getOneParamSet(p, 1)['w.init'][[1]]), "numeric")
  expect_equal(getAllParamSets(p)[, 'w.init'], c(4,5,6))
})

# Test getOneParamSet()
p = parameterSet()
my.params = getOneParamSet(p, 1)

test_that("getOneParamSet returns the correct number of values", {
  expect_equal(length(my.params), 14)
})

test_that("getOneParamSet return values are the proper class", {
  expect_equal(class(my.params$n.comp), "numeric")
  expect_equal(class(my.params$center.cols), "logical")
  expect_equal(class(my.params$scale.cols), "logical")
  expect_equal(class(my.params$center.rows), "logical")
  expect_equal(class(my.params$scale.rows), "logical")
  expect_equal(class(my.params$simultaneous), "logical")
  expect_equal(class(my.params$alg.typ), "character")
  expect_equal(class(my.params$fun), "character")
  expect_equal(class(my.params$alpha), "numeric")
  expect_equal(class(my.params$maxit), "numeric")
  expect_equal(class(my.params$tol), "numeric")
  expect_equal(class(my.params$w.init), "numeric")
  expect_equal(class(my.params$verbose), "logical")
})

