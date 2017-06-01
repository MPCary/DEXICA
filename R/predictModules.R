#' Predict gene modules using ICA
#'
#' @description Predicts gene modules (sets of co-expressed genes) from a previously
#' preprocessed (see Details) gene expression data matrix using independent component
#' analysis (ICA).
#'
#' @usage predictModules(x, n.comp = ceiling(ncol(x) / 4), ...)
#'
#' @param x a numeric matrix(-like object), typically with genes in rows and samples
#' (e.g., microarrays) in columns
#' @param n.comp an integer specifying the number of modules to predict
#' @param ... additional arguments to be passed to \code{\link[DEXICA]{dexFastICA}})
#'
#' @details This function runs an ICA algorithm (currently only fastICA is supported)
#' on the input matrix.  It does not conduct preprocessing operations (such as centering and
#' scaling on the rows and/or columns) on the input matrix first; \code{\link{preprocessMatrix}} should
#' be run on the input matrix to carryout such steps prior to using this function.
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{S}{The estimated source (gene module definition) matrix}
#'  \item{A}{The estimated mixing (e.g., weight of each module in each array) matrix}
#' }
#'
#' @seealso \code{\link[DEXICA]{dexFastICA}},
#' \code{\link[DEXICA]{preprocessMatrix}}
#'
#' @examples
#'x = matrix(rnorm(100), 10, 10)
#'x = preprocessMatrix(x)
#'m = predictModules(x, n.comp = 3)
#'
#' @include dexFastICA.R
#' @include checkIO.R
#'
#' @export

predictModules = function(x, n.comp = ceiling(ncol(x) / 4), ...) {
  res = dexFastICA(x, n.comp, ...)
  rownames(res$S) = rownames(x)
  colnames(res$A) = colnames(x)
  return(list(S = res$S, A = res$A, attempts = res$attempts, iterations = res$iterations,
              actual.tol = res$actual.tol, min.obs.tol = res$min.obs.tol,
              converged = res$converged))
}
