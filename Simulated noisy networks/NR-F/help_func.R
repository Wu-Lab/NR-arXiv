
#---------------------------------------------------------------------------------
#                                 help functions
#---------------------------------------------------------------------------------
# plot the membership of a graph
plot.membership <- function(graph, membership, main = "") {
  V(graph)$member <- membership
  mem.col <- c("#00CD00", "#FFB90F", '#00BFFF', '#EE799F')
  V(graph)$color <- mem.col[membership]
  plot(
    graph,
    edge.width = 2,
    vertex.color = V(graph)$color,
    main = main,
    edge.arrow.size = 0.7,
    vertex.label.cex = 0.8,
    vertex.size = 6
  )
}

# change the weighted matrix into the adjacency matrix
change_into_adj <- function(w, thre) {
  adj_mat <- matrix(0, nrow = n, ncol = n)
  adj_mat[w > thre] <- 1
  return(adj_mat)
}

# find the threshold such that the unweighted graph has no isolated vertices
find_thre <- function(W) {
  thre <- unique(sort(W, decreasing = T))
  for (k in n:length(thre)) {
    adj <- change_into_adj(W, thre[k])
    adj <- (adj + t(adj)) / 2
    rs <- rowSums(adj)
    if (length(rs[rs == 0]) == 0) {
      break
    }
  }
  
  return(adj)
}

# NR_F algorithm
NR_F <- function(A,m) {
  p1 <- A / rowSums(A)
  p2 <- (m-1)*p1 %*% solve(diag(m, nrow = n) - p1) # row stochastic matrix
  stationary_d <- Null(p2 - diag(1, nrow = n))
  stationary_d <- stationary_d / sum(stationary_d)
  if (dim(stationary_d)[2] != 1) {
    w_out <- matrix(0, nrow = n, ncol = n)
    print('wrong')
  }
  else {
    w_out <- diag(as.vector(stationary_d)) %*% p2
  }
  
  
  return(w_out)
}


# generate the planted l-partition benchmark graphs randomly
generate_plp <- function(n1, l, p_in, p_out) {
  n <- n1 * l
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      A[i, j] <- sample(c(0, 1), 1, prob = c(1 - p_out, p_out))
    }
  }
  
  n2 <- 1
  for (m in 1:l) {
    for (i in n2:(m * n1)) {
      if (i == n) {
        break
      }
      for (j in (i + 1):(m * n1)) {
        if (i != j) {
          A[i, j] <- sample(c(0, 1), 1, prob = c(1 - p_in, p_in))
        }
      }
    }
    n2 <- n2 + n1
  }
  
  A <- (A + t(A)) / 2
  A[A > 0] <- 1
  A <- A - diag(diag(A))
}



# compute noisy proportion
compute_error <- function(n1, l, Adj) {
  n <- n1 * l
  ref_mode <- matrix(1,nrow=n1,ncol=n1)
  ref <- 1-as.matrix(bdiag(ref_mode,ref_mode,ref_mode,ref_mode,ref_mode))
  
  outer <- Adj * ref
  outer_edges <- length(outer[outer == 1]) / 2
  sum_edges <- sum(Adj) / 2
  inner_edges <- sum_edges - outer_edges
  
  return(outer_edges / inner_edges)
}
