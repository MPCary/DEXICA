.checkMatrix <- function(x) {
  # This function checks to make sure x is a properly formed numeric matrix or, if not,
  # coerces it into one if possible.
  # This function is not exported.

  d = dim(x)
  if(is.null(d) | (length(d) != 2)) {
    stop("x must be a 2D matrix")
  }

  check.for.factors = sapply(x, is.factor)
  if(any(check.for.factors)){
    stop("x must not contain factor data")
  }

  # Preserve names
  rn = rownames(x)
  cn = colnames(x)

  # Convert to numeric (in case of character data)
  y = apply(x, 2, as.numeric)
  dim(y) = d

  if(sum(is.na(y)) > 0) {
    stop("x must not contain missing values or values coerced to NA")
  }

  rownames(y) = rn
  colnames(y) = cn
  return(y)
}

.check.xa = function(x, a) {
  x2 = .checkMatrix(x)
  a2 = .checkMatrix(a)
  common.genes = intersect(rownames(x2), rownames(a2))
  lcg = length(common.genes)
  if(lcg == 0) stop("Found no common IDs between x and a.")
  smaller = min(c(nrow(x), nrow(a)))
  if(lcg < (0.05 * smaller)) warning("x and a have less than 5% of IDs in common.")
  return(0)
}
