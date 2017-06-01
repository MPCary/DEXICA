#' Evaluate a set of modules using gene annotations
#'
#' @description Evaluates a set of gene modules (e.g., the output from
#' \code{\link{predictModules}}) by checking them for enrichment of gene
#' annotations.
#'
#' @usage checkModules(x, a, p = 0.05, return.p.matrix = FALSE)
#'
#' @param x a binary matrix indicating module membership (genes in rows, modules in columns)
#' @param a a binary gene annotation matrix (genes in rows, annotations in columns)
#' @param p FDR-adjusted p-value below which to consider an annotation significant
#' @param return.p.matrix logical value indicating whether the full p-value matrix should
#' be returned in the output
#'
#' @details This function evaluates a set of predicted gene modules by testing whether
#' each module is enriched for any gene annotations.  Enrichment p-values are calculated using
#' hypergeometric statistics (see \code{\link[stats]{phyper}}).  If x is not a binary matrix,
#' module membership is determined using \code{\link[DEXICA]{partitionComponents}} with default
#' paramaters.  The FDR-correction (see \code{\link[stats]{p.adjust}}) is applied to all
#' p-values prior to determining signficance.
#'
#' @return A list with the following elements:
#' \describe{
##'  \item{num.annotations}{Number of annotations tested (number of columns in a)}
##'  \item{num.modules}{Number of modules tested (number of columns in x)}
##'  \item{anns.signif}{Number of annotations that were significant in at least
##'  one module}
##'  \item{mods.signif}{Number of modules that were significant for at least
##'  one annotation}
##'  \item{p.matrix}{The module-by-annotation p-value matrix (if return.p.matrix == TRUE)}
##' }
#'
#' @examples
#'x = matrix(rnorm(100), 10, 10)
#'a = matrix(sample(c(0,1), size = 100, replace = TRUE), 10, 10)
#'m = predictModules(x)$S
#'rownames(a) = rownames(m) = as.character(1:10)
#'m.p = partition(m, t = 2)
#'m.p2 = splitMatrix(m.p)
#'checkModules(m.p2, a)
#'
#' @include checkIO.R
#' @export

checkModules = function(x, a, p = 0.05, return.p.matrix = FALSE) {
  # Check x and a
  check = .check.xa(x, a)
  if(check != 0) stop("x and a must be matrices with similar row names.")

  # Check if x is binary
  if(all(x %in% c(0,1))) {
    # Calculate phyper statistics
    result = .get.phyper.matrix(x, a)
  } else {
    # Paritition and split x
    x.p = partition(x)
    x.p.split = splitMatrix(x.p)
    result = .get.phyper.matrix(x.p.split, a)
  }

  # Adjust p-values
  result[] = p.adjust(result, method = "fdr")

  # Count significant annotations / modules
  # Turning off warnings because usually some modules have no genes (after partitioning)
  # Note: Such modules do not impact p-value adjustments due to their p-values being NA
  oldw <- getOption("warn")
  options(warn = -1)

  signif.anns = sum(apply(result, 1, FUN = function(x) min(x, na.rm = TRUE) <= p))
  signif.mods = sum(apply(result, 2, FUN = function(x) min(x, na.rm = TRUE) <= p))

  # Turn warnings back on
  options(warn = oldw)

  # Generate return values
  if(return.p.matrix) {
    return(list(num.annotations = nrow(result),
                num.modules = ncol(result),
                anns.signif = signif.anns,
                mods.signif = signif.mods,
                p.matrix = result))
  } else {
    return(list(num.annotations = nrow(result),
                num.modules = ncol(result),
                anns.signif = signif.anns,
                mods.signif = signif.mods))
  }
}

.get.phyper.matrix = function(x, a) {
  # Calculate p-value matrix from x and a

  # Reduce matrices to common rows
  common.rows = intersect(rownames(x), rownames(a))
  x2 = x[common.rows, ]
  a2 = a[common.rows, ]

  # Create cross matrix
  cross = t(a2) %*% x2
  cross.dim = dim(cross)

  # Get hypergeometric test parameter k
  k = apply(x2, 2, sum) # Number of genes per module
  k.mat = matrix(k, ncol = cross.dim[2], nrow = cross.dim[1], byrow = TRUE)
  zero.cols = k == 0

  # Get hypergeometrix test parameters m and n
  m = apply(a2, 2, sum) # Number of genes per annotation
  m.mat = matrix(m, ncol = cross.dim[2], nrow = cross.dim[1], byrow = FALSE)
  n = dim(a2)[1] - m # Number of genes not annotated with each annotation
  n.mat = matrix(n, ncol = cross.dim[2], nrow = cross.dim[1], byrow = FALSE)

  # phyper() will give the lower tail, which is p(X) > x; we need p(X) >= x,
  # which we can get if we subtract 1 first
  cross.adj = cross - 1

  # Perform hypergeometric test
  phyper.mat = phyper(cross.adj, m.mat, n.mat, k.mat, lower.tail = FALSE)

  # No tests were performed on the zero columns; set these to NA
  phyper.mat[, zero.cols] = NA

  # Set names and return the p-value matrix
  colnames(phyper.mat) = colnames(x2)
  rownames(phyper.mat) = colnames(a2)
  return(phyper.mat)
}

# .get.ks.matrix = function(x, a) {
#   stop("Sorry, KS test p-values have not been implemented yet (check future versions of this package).")
# }
