#' Split matrix
#'
#' Given a numeric matrix, this function performs a column-wise separation of positive
#' values from negative values, producing a return matrix with twice the number of columns
#' as the input matrix.
#'
#' @usage splitMatrix(x)
#'
#' @param x A numeric matrix
#'
#' @details Conceptually, this function makes two copies of the input matrix (we will call
#' them A and B).  It replaces all negative values in copy A with 0, and all positive
#' values in copy B with 0.  It then multiplies copy B by -1, and finally combines (via \code{\link{cbind}})
#' the two copies together column-wise.  The strings "_pos" and "_neg" are appended to the
#' original column names to generate the column names for the new matrix.
#'
#' @return A matrix devoid of negative values with twice the number of columns as
#' the input matrix.
#'
#' @examples
#' x = matrix(rnorm(9), 3, 3)
#' splitMatrix(x)
splitMatrix <- function(x) {
  x = .checkMatrix(x) # Check x for well-formedness

  # Sort out the names
  n = ncol(x)
  old_col_names = colnames(x)
  if(is.null(old_col_names)) old_col_names = as.character(1:n)
  pos_names = paste(old_col_names, "_pos", sep = "")
  neg_names = paste(old_col_names, "_neg", sep = "")
  new_names = c(pos_names, neg_names)

  # Copy x, then alter it
  b = x * -1
  b[b < 0] = 0
  x[x < 0] = 0

  # Create return matrix
  b2 = cbind(x, b)
  colnames(b2) = new_names
  rownames(b2) = rownames(x)

  return(b2)
}
