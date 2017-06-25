#' Preprocess a data matrix
#'
#' @description Centers and/or scales the columns and/or rows of the input matrix.
#'
#' @usage preprocessMatrix(x, center.cols = TRUE, scale.cols = FALSE, center.rows = FALSE,
#' scale.rows = FALSE, simultaneous = FALSE, verbose = TRUE)
#'
#' @param x a numeric matrix(-like object), typically with genes in rows and samples
#' (e.g., microarrays) in columns
#' @param center.cols a logical value indicating whether to adjust columns to have 0 mean
#' @param scale.cols a logical value indicating whether to adjust columns to have unit variance
#' @param center.rows a logical value indicating whether to adjust rows to have 0 mean
#' @param scale.rows a logical value indicating whether to adjust rows to have unit variance
#' @param simlutanous a logical value indicating whether to adjust rows and columns simultaneously
#' @param verbose a logical value indicating the desired level of output as the algorithm runs
#'
#' @details This function centers and/or scales the columns and/or rows of the input matrix.
#' Default parameter values correspond to those of the fastICA algorithm.
#'
#' @return A centered and/or scaled matrix.
#'
#' @seealso \code{\link[base]{scale}},
#' \code{\link[scaleMatrix]{scaleMatrix}}
#'
#' @examples
#'x = matrix(rnorm(100), 10, 10)
#'x = preprocessMatrix(x)
#'
#' @export

preprocessMatrix <- function(x, center.cols = TRUE, scale.cols = FALSE, center.rows = FALSE,
                             scale.rows = FALSE, simultaneous = FALSE, verbose = FALSE) {
  # Center and/or scale
  if(!simultaneous) {
    if(verbose) message("Centering and scaling...")
    x = scale(x, center = center.cols, scale = scale.cols)
    x = t(scale(t(x), center = center.rows, scale = scale.rows))
  } else {
    x = scaleMatrix(x, center.cols = center.cols, scale.cols = scale.cols,
                    center.rows = center.rows, scale.rows = scale.rows,
                    verbose = verbose)
  }
  return(x)
}
