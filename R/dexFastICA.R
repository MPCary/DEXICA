#' Run the fastICA algorithm without matrix preprocessing
#'
#' @description Runs the fastICA algorithm but does not perform any matrix preprocessing operations
#' first (this differs from the \code{fastICA} function in the \pkg{fastICA} package.)  Preprocessing steps should be carried out
#' by running other functions, e.g. \code{\link{preprocessMatrix}}, on the input matrix prior to
#' using this function.
#'
#' @usage dexfastICA(X, n.comp, alg.typ = c("parallel", "deflation"), fun = c("logcosh", "exp"),
#' alpha = 1, maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL, max.attempts = 3)
#'
#' @param X a numeric matrix(-like object)
#' @param n.comp an integer specifying the number of components to extract
#' @param alg.typ,fun,alpha,maxit,tol,verbose See \code{\link[fastICA]{fastICA}} documentation
#' @param w.init integer value to reset random number generator with (see Details)
#' @param max.attempts maximum number of attempts at convergence to make
#'
#' @details This function runs the fastICA algorithm on the input matrix but it does so without running
#' any preprocessing operations (such as centering and scaling on the rows and/or columns) on the input
#' matrix first.  This allows alternative preprocessing methods (such as \code{\link[scaleMatrix]{scaleMatrix}})
#' to be applied to the input matrix, which can improve results.  The \code{\link{preprocessMatrix}} function
#' may be run on the input matrix to carryout typical preprocessing steps, such as column centering.
#'
#' Note: The default parameter values in \code{\link{preprocessMatrix}} match those of \code{\link[fastICA]{fastICA}},
#' thus, running \code{preprocessMatrix} followed by \code{dex.fastICA} with default parameter values is equivalent to running \code{fastICA}.
#'
#' If \code{w.init} is not NULL, its value will be used to reset the random number generator
#' (using \code{\link[base]{set.seed}}) prior to randomly generarting the initial W ('unmixing') matrix.
#' This differs from \code{\link[fastICA]{fastICA}}, as w.init there can be either NULL or a matrix.
#' This modification to w.init was made in order to make reproducing results easier.
#'
#' Another difference between this algorithm and \code{\link[fastICA]{fastICA}} is that multiple attempts
#' at convergence can be made if convergence fails in the first attempt.  The maximum number of attempts
#' that will be made is specified by the \code{max.attempts} parameter.  The maximum number of iterations in
#' each attempt, specified by \code{maxit}, is increased by 1.5x in each subsequent attempt.  If convergence
#' fails, the ICA solution found at the last iteration of the last attempt will be returned as an approximate
#' solution.
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{S}{The estimated source matrix}
#'  \item{A}{The estimated mixing matrix}
#'  \item{attempts}{The number of attempts made at finding a solution}
#'  \item{iterations}{The number of iterations made during the last attempt}
#'  \item{converged}{A logical value indicating whether convergence was reached}
#' }
#'
#' @seealso \code{\link[fastICA]{fastICA}}, \code{\link[DEXICA]{preprocessMatrix}},
#' \code{\link[DEXICA]{predictModules}}
#'
#' @examples
#'x = matrix(rnorm(100), 10, 10)
#'x = preprocessMatrix(x)
#'m = dexFastICA(x, n.comp = 3)
#'m$converged
#'
#'@include dex.ica.R.def.R
#'@include dex.ica.R.par.R
#'
#' @export
dexFastICA = function(X, n.comp, alg.typ = c("parallel", "deflation"),
                      fun = c("logcosh", "exp"), alpha = 1,
                      maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL,
                      max.attempts = 3) {
  # Check input data matrix
  X = DEXICA:::.checkMatrix(X)

  # Check parameters
  if (alpha < 1 || alpha > 2)
    stop("alpha must be in range [1,2]")
  alg.typ.match <- match.arg(alg.typ)
  fun.match <- match.arg(fun)
  n <- nrow(X)
  p <- ncol(X)
  if (n.comp > min(n, p)) {
    message("'n.comp' is too large: reset to ", min(n, p))
    n.comp <- min(n, p)
  }
  if (is.null(w.init))
    w.init.matrix <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  else {
    set.seed(w.init)
    w.init.matrix <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  }

  # Note: Only the R version of fastICA is supported
  # Whiten data matrix
  X = t(X)
  if (verbose) message("Whitening...")
  V <- X %*% t(X)/n
  s <- La.svd(V)
  D <- diag(c(1/sqrt(s$d)))
  K <- D %*% t(s$u)
  K <- matrix(K[1:n.comp, ], n.comp, p)
  X1 <- K %*% X

  # Attempt ICA; repeat if convergence fails
  attempt = 1
  while(attempt <= max.attempts) {
    if(verbose) message("Running ICA attempt ", attempt, "...")
    attempt = attempt + 1
    tryCatch({
      my.run = if (alg.typ.match == "deflation") {
        dex.ica.R.def(X1, n.comp, tol = tol, fun = fun.match, alpha = alpha,
                  maxit = maxit + 1, verbose = verbose, w.init = w.init.matrix)
      } else if (alg.typ.match == "parallel") {
        dex.ica.R.par(X1, n.comp, tol = tol, fun = fun.match, alpha = alpha,
                  maxit = maxit + 1, verbose = verbose, w.init = w.init.matrix)
      } else stop("Unknown alg.typ")
    }, error = function(e) {
      message("ICA_RUN_ERROR: ", e)
    })

    # Get results
    a = my.run$W
    last.it = my.run$iterations - 1 # Last iteration
    actual.tol = my.run$actual.tol
    min.obs.tol = my.run$min.obs.tol

    # Decide whether to try again
    if(actual.tol <= tol) {
      if(verbose) message("Convergence achieved in ", last.it, " iterations.")
      w <- a %*% K
      S <- w %*% X
      A <- t(w) %*% solve(w %*% t(w))
      result = list(
        A = t(A),
        S = t(S),
        attempts = attempt - 1,
        iterations = last.it,
        actual.tol = actual.tol,
        min.obs.tol = min.obs.tol,
        converged = TRUE)
      return(result)
    } else {
      # Convergence failed, try again
      if(verbose) message("Failed to converge after ", last.it, " iterations.")
      maxit = round(maxit * 1.5)
    }
  }
  if(verbose) message("Convergence failed, returning best estimate.")
  w <- a %*% K
  S <- w %*% X
  A <- t(w) %*% solve(w %*% t(w))
  result = list(
    A = t(A),
    S = t(S),
    attempts = attempt - 1,
    iterations = last.it,
    actual.tol = actual.tol,
    min.obs.tol = min.obs.tol,
    converged = FALSE)
  return(result)
}
