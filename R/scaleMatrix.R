#' Scale rows and columns of a matrix
#'
#' @description Centers and/or scales the rows and/or columns of a numeric matrix.
#'
#' @usage scaleMatrix(x, center.rows = TRUE, center.cols = TRUE, scale.rows = TRUE,
#' scale.cols = TRUE, steps = 1000, bad.row.action = c("resample", "remove"),
#' bad.col.action = c("resample", "remove"), set.na.to = 0, verbose = FALSE)
#'
#' @param x a numeric matrix(-like object)
#' @param center.rows a logical value indicating whether to set row means to 0
#' @param center.cols a logical value indicating whether to set column means to 0
#' @param scale.rows a logical value indicating whether to set row variances to 1
#' @param scale.cols a logical value indicating whether to set column variances to 1
#' @param steps an integer
#' @param bad.row.action action to take when a row in x has 0 variance
#' @param bad.col.action action to take when a column in x has 0 variance
#' @param set.na.to numeric value to set NA values in x to
#' @param verbose a logical value indicating whether to show status messages during execution
#'
#' @details This function is a generalization of the generic function \code{\link[base]{scale}}, which
#' operates on only the columns of a matrix.  This function (optionally) operates on both
#' the rows and columns simultaneously, providing a resulting matrix with both rows and
#' columns having zero means and unit variance.
#'
#' The algorithm converges toward the solution in a stepwise manner - the larger the number
#' of steps, the more accurate (i.e., the closer to the desired outcome) will be the final
#' solution.  A step size of 1000 resulted in convergence (deviation from target < 2.2e-16)
#' in most test cases, but may be time consuming on large matrices.  Lower values for
#' \code{steps} will execute faster, but with possible loss of accuracy; a \code{step} size
#' < 100 is not recommended.
#'
#' If either: a) both \code{center.cols} and \code{scale.cols} are \code{FALSE}, or b) both
#' \code{center.rows} and \code{scale.rows} are \code{FALSE}, the function defaults to the
#' generic \code{\link[base]{scale}} function (using transposition to handle the latter case).
#'
#' @return A scaled and/or centered matrix.
#'
#'@examples
#'x = matrix(rnorm(100), 10, 10)
#'range(apply(x, 1, mean)) # row means
#'range(apply(x, 1, sd)) # row sd
#'range(apply(x, 2, mean)) # column means
#'range(apply(x, 2, sd)) # column sds
#'
#'x.scaled = scaleMatrix(x)
#'range(apply(x.scaled, 1, mean)) # row means
#'range(apply(x.scaled, 1, sd)) # row sds
#'range(apply(x.scaled, 2, mean)) # column means
#'range(apply(x.scaled, 2, sd)) # column sds
#'
#' @export
scaleMatrix <- function(x, center.rows = TRUE, center.cols = TRUE,
                        scale.rows = TRUE, scale.cols = TRUE,
                        steps = 1000,
                        bad.row.action = c("resample", "remove"),
                        bad.col.action = c("resample", "remove"),
                        set.na.to = 0, verbose = FALSE) {

  # If only only rows (or only columns) are to be modified, use the default
  # function
  if(!center.rows & !scale.rows) {
    return(scale(x, center = center.cols, scale = scale.cols))
  }

  if(!center.cols & !scale.cols) {
    return(t(scale(t(x), center = center.rows, scale = scale.rows)))
  }

  # NA values need to be removed
  if(is.na(set.na.to) | !is.numeric(set.na.to)) stop("set.na.to must be a numeric value.")
  x[is.na(x)] = set.na.to

  # If scaling rows, find bad rows (rows with zero variance)
  if(scale.rows) {
    x = .scale.rows(x, bad.row.action[1], verbose)
  }

  # If scaling cols, find bad cols (cols with zero variance)
  if(scale.cols) {
    x = .scale.cols(x, bad.col.action[1], verbose)
  }

  for(i in c(1:steps)) {
    x.r = t(scale(t(x), center.rows, scale.rows))
    x.c = scale(x, center.cols, scale.cols)
    x.r.d = x - x.r
    x.c.d = x - x.c
    mean.diff = (x.r.d + x.c.d) / 2
    mean.diff = mean.diff * (i/steps)
    if(verbose) message("Iteration ", i, " mean / max deviation from target: ", mean(mean.diff),
                        " / ", max(abs(mean.diff)))
    if(max(abs(mean.diff)) < 2.2e-16) {
      if(verbose) message("Convergence reached at iteration ", i)
      break
    }
    x = x - mean.diff
  }
  y = matrix(as.vector(x), dim(x)) # strip attributes
  dimnames(y) = dimnames(x)
  #return(matrix(as.vector(x), dim(x)))
  return(y)
}

.scale.rows = function(x, bad.row.action, verbose) {
  if(verbose) message("Checking rows... ", appendLF = FALSE)
  var.r = apply(x, 1, var, na.rm = TRUE)
  row.check = sum(var.r == 0)
  if(row.check == 0 & verbose) message("all rows OK.") else {
    if(verbose) {
      message("found ", row.check, " rows with zero variance, will ", bad.row.action, " them.")
    }
    x = switch(bad.row.action,
               remove = x[var.r != 0, ],
               resample = .resample.rows(x, var.r),
               stop("Unrecognized bad.row.action"))
  }
  return(x)
}

.resample.rows = function(x, var.r) {
  # Resamples bad rows from the rest of the matrix
  x[var.r == 0, ] = sample(x[var.r != 0], size = (ncol(x) * sum(var.r == 0)),
                         replace = TRUE)
  return(x)
}

.scale.cols = function(x, bad.col.action, verbose) {
  if(verbose) message("Checking cols... ", appendLF = FALSE)
  var.c = apply(x, 2, var, na.rm = TRUE)
  col.check = sum(var.c == 0)
  if(col.check == 0 & verbose) message("all columns OK.") else {
    if(verbose) {
      message("found ", col.check, " columns with zero variance, will ", bad.col.action[1], " them.")
    }
    x = switch(bad.col.action,
               remove = x[, var.c != 0],
               resample = .resample.cols(x, var.c),
               stop("Unrecognized bad.col.action"))
  }
  return(x)
}

.resample.cols = function(x, var.c) {
  # Resamples bad cols from the rest of the matrix
  x[, var.c == 0] = sample(x[, var.c != 0], size = (nrow(x) * sum(var.c == 0)),
                           replace = TRUE)
  return(x)
}
