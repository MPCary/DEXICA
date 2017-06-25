library(DEXICA)
context("Test batch mode functionality")

# Prepare test compendia
set.seed(1)
x = matrix(rnorm(100), 10, 10)
y = matrix(rnorm(100), 10, 10)

# Prepare annotation matrices
a1 = matrix(sample(c(0,1), size = 100, replace = TRUE), 10, 10)
a2 = matrix(sample(c(0,1), size = 100, replace = TRUE), 10, 10)
a3 = matrix(sample(c(0,1), size = 100, replace = TRUE), 10, 10)

# Set rownames
rownames(x) = rownames(y) = rownames(a1) = rownames(a2) = rownames(a3) = as.character(1:10)

# Create dexBatch object
db = dexBatch(compen = list("x" = x, "y" = y), annmats = list("a1" = a1, "a2" = a2, "a3" = a3))

test_that("the dexBatch object was created successfully", {
  expect_equal(class(db)[[1]], "DexBatch")
})

test_that("countJobs() works properly", {
  expect_equal(countJobs(db), 2)
})

# Create with a custom ParameterSet
p = parameterSet(w.init = c(1,2,3))
db = dexBatch(compen = list("x" = x, "y" = y), annmats = list("a1" = a1, "a2" = a2, "a3" = a3),
              params = p)

test_that("the dexBatch object was created successfully", {
  expect_equal(class(db)[[1]], "DexBatch")
  expect_equal(countJobs(db), 6)
})

# Test compatibility with LazyData (if that's installed)
if("LazyData" %in% installed.packages()) {
  library(LazyData)
  p = parameterSet(n.comp = c(5,10,20))
  x = lazyMatrix(name = "GPL200.50.mini", package = "DEXICA")
  db = dexBatch(compen = list("x" = x, "y" = y), annmats = list("a1" = a1, "a2" = a2, "a3" = a3),
                params = p)

  test_that("the dexBatch object was created successfully", {
    expect_equal(class(db)[[1]], "DexBatch")
    expect_equal(countJobs(db), 6)
  })
}

# Test connection
temp = tempfile()
db = dexBatch(compen = list("x" = x, "y" = y), annmats = list("a1" = a1, "a2" = a2, "a3" = a3),
              params = p, output = "test.results")

test_that("the dexBatch object can write to a file", {
  con = file(db@output)
  expect_true("connection" %in% class(con))
  close(con)
})

# Test runJob
# db = dexBatch(compen = list("x" = x, "y" = y), annmats = list("a1" = a1, "a2" = a2, "a3" = a3),
#               params = p)
# runJob(db, 1)
