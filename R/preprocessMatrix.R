#' Preprocess a data matrix
#'
#' @description Centers and/or scales the columns and/or rows of the input matrix.
#'
#' @usage preprocessMatrix(x, center.cols = TRUE, scale.cols = FALSE, center.rows = FALSE,
#' scale.rows = FALSE, simultaneous = FALSE, set.na.to = c("mean", "median"), verbose = TRUE)
#'
#' @param x a numeric matrix(-like object), typically with genes in rows and samples
#' (e.g., microarrays) in columns
#' @param center.cols a logical value indicating whether to adjust columns to have 0 mean
#' @param scale.cols a logical value indicating whether to adjust columns to have unit variance
#' @param center.rows a logical value indicating whether to adjust rows to have 0 mean
#' @param scale.rows a logical value indicating whether to adjust rows to have unit variance
#' @param simlutanous a logical value indicating whether to adjust rows and columns simultaneously
#' @param set.na.to character string indicating how to handle `NA` values
#'
#' @details This function centers and/or scales the columns and/or rows of the input matrix.
#' Default parameter values correspond to those of the fastICA algorithm in the fastICA
#' package.
#'
#' If the matrix contains `NA` values after centering/scaling (this can happen, e.g., if a row
#' in `x` has 0 variance), the `NA` values are replaced using the value indicated by `set.na.to`
#' (either the global mean or median.)
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
                             scale.rows = FALSE, simultaneous = FALSE,
                             set.na.to = c("resample", "mean", "median"), verbose = FALSE) {
  # Check set.na.to
  set.na.to = set.na.to[[1]]
  if(!(set.na.to %in% c("mean", "median", "resample"))) stop("Unrecognized set.na.to value.")

  # Center and/or scale
  if(!simultaneous) {
    if(verbose) message("Centering and scaling...")
    x = scale(x, center = center.cols, scale = scale.cols)
    x = t(scale(t(x), center = center.rows, scale = scale.rows))
    x[is.na(x)] = switch(set.na.to,
                         mean = mean(x, na.rm = TRUE),
                         median = median(x, na.rm = TRUE),
                         resample = {
                           if(length(x[!is.na(x)]) == 0) stop("No non-NA values to sample from")
                           sample(x[!(is.na(x))], sum(is.na(x)), replace = TRUE)
                         },
                         stop("Could not perform set.na.to action"))
  } else {
    x = scaleMatrix(x, center.cols = center.cols, scale.cols = scale.cols,
                    center.rows = center.rows, scale.rows = scale.rows,
                    set.na.to =  set.na.to, zero.var.action = "resample",
                    verbose = verbose)
  }
  return(x)
}
