
#---------------------------------------------------------------------------------
#                                 help functions
#---------------------------------------------------------------------------------

# Transform the weighted matrix "w" into adjacency matrix based on the threshold
change_into_adj <- function(w, thre) {
  adj_mat <- matrix(0, nrow = n, ncol = n)
  adj_mat[w > thre] <- 1
  return(adj_mat)
}


# NR_F algorithm
NR_F <- function(A) {
  p1 <- A / rowSums(A)
  p2 <- p1 %*% solve(diag(2, nrow = n) - p1)
  stationary_d <- Null(p2 - diag(1, nrow = n))
  stationary_d <- stationary_d / sum(stationary_d)
  w_out <- diag(as.vector(stationary_d)) %*% p2
  return(w_out)
}
