#---------------------------------------------------------------------------------
#                                 help functions
#---------------------------------------------------------------------------------

NR_B <- function(mat, m, eps1, eps2) {
  
  n <- dim(mat)[1]
  mat <- mat + eps1 + eps2 * diag(n)
  
  P1 <- mat / rowSums(mat)
  P2 <- m * solve((m - 1) * diag(n) + P1) %*% P1
  P2 <- P2 - pmin(apply(P2, 1, min), 0)
  P2 <- P2 / rowSums(P2)
  
  stat_d <- abs(Null(P2 - diag(n)))
  net_new <- diag(as.vector(stat_d)) %*% P2
  
  net_new <- net_new + t(net_new)
  output_network <- net_new - diag(diag(net_new))
  output_network <-
    (output_network - min(output_network)) / (max(output_network) - min(output_network))
  
}
