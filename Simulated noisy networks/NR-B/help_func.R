


change_into_adj_keep_edges <- function(w,edges){
  temp <- sort(w,decreasing = T)
  thre <- temp[2*edges]
  A <- matrix(0,nrow = n,ncol = n)
  A[w>=thre] <- 1
  A <- A+t(A)
  A[A>0] <- 1
  return(A)
}



change_into_adj <- function(w,thre){
  adj_mat <- matrix(0,nrow = n,ncol = n)
  adj_mat[w>thre] <- 1
  return(adj_mat)
}


pre_processing <- function(A,m=2){
  n <- dim(A)[1]
  A1 <- A
  for (eps_2 in seq(1,5,0.5)){
    for (eps_1 in seq(1,5,0.5)){
      A <- A1+eps_1*matrix(1,nrow = n,ncol = n)+eps_2*diag(n)
      P1 <- A/rowSums(A)
      P2 <- solve((m-1)*diag(1,nrow = n)+P1)%*%(m*P1)
      if (length(P2[P2<0])==0){
        print('break')
        result <- list(P2,eps_1,eps_2)
        return(result)
      }
    }
  }
}

# compute noisy proportion
compute_error <- function(A_true,A_pre) {
  error_edges <- sum((A_true==0)&(A_pre==1))
  true_edges <- sum((A_true==1)&(A_pre==1))
    
  SNR_inverse <- error_edges/true_edges
  return(SNR_inverse)

}
