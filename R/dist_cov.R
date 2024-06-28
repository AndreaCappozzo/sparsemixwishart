#' Compute Non-Euclidean distances for covariance matrices, according to Dryden AOAS 2009
#' Useful for subsequently performing hierarchical clustering with different distance types
#'
#' @param data
#' @param method
#'
#' @return
#' @export
#'
#' @examples
dist_cov_creator <- function(data, method = "Euclidean") {
  N <- dim(data)[3]
  dist_matrix <- matrix(NA, nrow = N, ncol = N)

  for (row in 1:N) {
    for (col in row:N) {
      dist_row_col <- shapes::distcov(
        S1 = data[, , row],
        S2 = data[, , col],
        method = method
      )
      dist_matrix[row, col] <- dist_row_col
      dist_matrix[col, row] <- dist_row_col
    }
  }
  stats::as.dist(dist_matrix)
}
