#---------------------------------------------------------------------------------
'''
Parameters Explanation:
mat: input noisy matrix
m: Parameters of NR_B, controlling the extent of matrix modification
eps1: Value added to the whole matrix in preprocessing step
eps2: Value added to the diagonal in preprocessing step
'''
#---------------------------------------------------------------------------------

NR_B <- function(mat, m, eps1, eps2) {
  
  # preprocessing
  n <- dim(mat)[1]
  mat <- mat + eps1 + eps2 * diag(n)
  
  # change into transition matrix
  P1 <- mat / rowSums(mat)
  P2 <- m * solve((m - 1) * diag(n) + P1) %*% P1
  
  # minus the smallest negative value from each row with negative numbers
  P2 <- P2 - pmin(apply(P2, 1, min), 0)
  P2 <- P2 / rowSums(P2)
  
  
  # compute the stationary distribution and the output matrix
  stat_d <- abs(Null(P2 - diag(n)))
  net_new <- diag(as.vector(stat_d)) %*% P2
  
  net_new <- net_new + t(net_new)
  output_network <- net_new - diag(diag(net_new))
  output_network <-
    (output_network - min(output_network)) / (max(output_network) - min(output_network))
  
}
