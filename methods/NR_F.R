#---------------------------------------------------------------------------------
'''
Parameters Explanation:
mat: input noisy matrix
m: Parameters of NR_F, controlling the extent of matrix modification
'''
#---------------------------------------------------------------------------------
library(MASS)

NR_F <- function(mat, m) {
  # change into transition matrix
  n <- dim(mat)[1]
  p1 <- mat / rowSums(mat)
  
  # diffusion
  p2 <- (m - 1) * p1 %*% solve(diag(m, nrow = n) - p1)
  
  # compute the stationary distribution and the output matrix
  stationary_d <- Null(p2 - diag(1, nrow = n))
  stationary_d <- stationary_d / sum(stationary_d)
  if (dim(stationary_d)[2] == 0) {
    print('The stationary distribution does not exist！')
  } else if (dim(stationary_d)[2] > 1) {
    print('The stationary distribution is not unique！')
  } else{
    w_out <- diag(as.vector(stationary_d)) %*% p2
    return(w_out)
  }
}
