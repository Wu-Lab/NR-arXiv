library(igraph)
library(ggplot2)
library(expm)
library(MASS)
source('help_func.R')

#-----------------------------------------------------------------------------------
# NR-F improves the accuracy of Community Detection of Zachary karate club network.
#-----------------------------------------------------------------------------------
# We evaluate the clustering performance of six community detection methods 
# (fast-greedy, louvain, optimal, leading-eigen, Walktrap, spinglass) 
# on the Zachary network before and after applying NR-F and NE as a preprocessing 
# step once and twice respectively.
#-----------------------------------------------------------------------------------


Methods <- c(
  'cluster_fast_greedy',
  'cluster_louvain',
  'cluster_optimal',
  'cluster_leading_eigen',
  'cluster_walktrap',
  'cluster_spinglass'
)

g_origi <- make_graph("Zachary")
true_membership <- c(1,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1,
                     1,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2)

n <- length(V(g_origi))
e1 <- length(E(g_origi))

adj_origi <- matrix(as_adjacency_matrix(g_origi), nrow = n)


#-------------------------------------------------------------------------------
# apply NR_F algorithm once
print('Applying NR-F once................')
out_NR_F <- NR_F(adj_origi)
out_NR_F <- (out_NR_F + t(out_NR_F)) / 2
out_NR_F <- out_NR_F - diag(diag(out_NR_F))
out_NR_F <- out_NR_F / max(out_NR_F)

adj_NR_F <- change_into_adj(out_NR_F, 0)
g_NR_F <-
  graph_from_adjacency_matrix(adj_NR_F, mode = 'undirected', weighted = T)
E(g_NR_F)$weight <- out_NR_F[lower.tri(out_NR_F)]

#-------------------------------------------------------------------------------
# apply NR_F algorithm twice
print('Applying NR-F twice................')
out_NR_F_twice <- NR_F(out_NR_F)
out_NR_F_twice <- (out_NR_F_twice + t(out_NR_F_twice)) / 2
out_NR_F_twice <- out_NR_F_twice - diag(diag(out_NR_F_twice))
out_NR_F_twice <- out_NR_F_twice / max(out_NR_F_twice)

adj_NR_F_twice <- change_into_adj(out_NR_F_twice, 0)
g_NR_F_twice <-
  graph_from_adjacency_matrix(adj_NR_F_twice, mode = 'undirected', weighted = T)
E(g_NR_F_twice)$weight <- out_NR_F_twice[lower.tri(out_NR_F_twice)]


#-------------------------------------------------------------------------------
# apply NE algorithm once
print('Applying NE once................')
out_NE_once <-
  as.matrix(read.csv(
    "zachary_NE.csv",
    sep = ',',
    encoding = "UTF-8",
    header = F
  ))

out_NE_once <- (out_NE_once + t(out_NE_once)) / 2
out_NE_once <- out_NE_once - diag(diag(out_NE_once))
out_NE_once <- out_NE_once / max(out_NE_once)

adj_NE_once <- change_into_adj(out_NE_once, 0)
g_NE_once <-
  graph_from_adjacency_matrix(adj_NE_once, mode = 'undirected', weighted = T)
E(g_NE_once)$weight <- out_NE_once[lower.tri(out_NE_once)]

#-------------------------------------------------------------------------------
# apply NE algorithm twice
print('Applying NE twice................')
out_NE_twice <-
  as.matrix(read.csv(
    "zachary_NE2.csv",
    sep = ',',
    encoding = "UTF-8",
    header = F
  ))

out_NE_twice <- (out_NE_twice + t(out_NE_twice)) / 2
out_NE_twice <- out_NE_twice - diag(diag(out_NE_twice))
out_NE_twice <- out_NE_twice / max(out_NE_twice)

adj_NE_twice <- change_into_adj(out_NE_twice, 0)
g_NE_twice <-
  graph_from_adjacency_matrix(adj_NE_twice, mode = 'undirected', weighted = T)
Wei <- out_NE_twice[lower.tri(out_NE_twice)]
E(g_NE_twice)$weight <- Wei[Wei>0]



print('*********************************************************************')
print('evaluating the clustering performance.........')
score_nmi <- matrix(0, nrow = length(Methods), ncol = 5)
score_rand <- matrix(0, nrow = length(Methods), ncol = 5)
score_ARI <- matrix(0, nrow = length(Methods), ncol = 5)
col_content <-
  c('g_origi', 'g_NE_once', 'g_NE_twice', 'g_NR_F', 'g_NR_F_twice')



for (i in 1:length(Methods)) {
  print(paste0('Methods: ',Methods[i]))
  if (i < 6) {
    for (j in 1:5) {
      Compare_membership <- get(Methods[i])(get(col_content[j]),
                                            weights = E(get(col_content[j]))$weight)$membership
      score_nmi[i, j] <-
        compare(true_membership, Compare_membership,
                method = 'nmi')
      
      score_ARI[i, j] <-
        compare(true_membership, Compare_membership,
                method = 'adjusted.rand')
      
      score_rand[i, j] <-
        compare(true_membership, Compare_membership,
                method = 'rand')
    }
    
  }
  else {
    #cluster_spinglass gets different result every time,
    # so we take the average of 20 trials
    for (loop in 1:20) {
      for (j in 1:5) {
        Compare_membership <- get(Methods[i])(get(col_content[j]),
                                              weights = E(get(col_content[j]))$weight)$membership
        score_nmi[i, j] <-
          score_nmi[i, j] + compare(true_membership, Compare_membership,
                                    method = 'nmi') / 20
        
        score_ARI[i, j] <-
          score_ARI[i, j] + compare(true_membership, Compare_membership,
                                    method = 'adjusted.rand') / 20
        
        score_rand[i, j] <-
          score_ARI[i, j] + compare(true_membership, Compare_membership,
                                    method = 'rand') / 20
        
      }
      
    }
  }
  
}
print('*********************************************************************')


#-------------------------------------------------------------------------------
#plot
methods_name <- c('greedy',
                  'louvain',
                  'optimal',
                  'eigen',
                  'walktrap',
                  'spinglass')

data <- data.frame(
  algorithm = rep(methods_name, 5),
  methods = c(
    rep('raw', length(methods_name)),
    rep('NE', length(methods_name)),
    rep('NE*2', length(methods_name)),
    rep('NR-F', length(methods_name)),
    rep('NR-F*2', length(methods_name))
  ),
  NMI = as.vector(score_nmi),
  ARI = as.vector(score_ARI),
  rand = as.vector(score_rand)
)


data$algorithm <- factor(data$algorithm,
                         c(
                           'greedy',
                           'louvain',
                           'optimal',
                           'eigen',
                           'walktrap',
                           'spinglass'
                         ))

data$methods <-
  factor(data$methods, c('raw', 'NE', 'NE*2', 'NR-F', 'NR-F*2'))




#-------------------------------------------------------------------------------
# plot result
ggplot(data, mapping = aes(x = algorithm, y = NMI, fill = methods)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#4682B4","#F0E442", "#E69F00")) +
  ylab('NMI') + xlab("") + theme_bw() + ylim(0, 1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,family='serif'),
        axis.text.x = element_text(size=13,family='serif',
                                   angle = 45,hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=13,family='serif'),
        axis.title.y = element_text(size=13,family='serif'))


ggplot(data, mapping = aes(x = algorithm, y = ARI, fill = methods)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#4682B4", "#F0E442", "#E69F00")) +
  ylab('ARI') + xlab("") + theme_bw() + ylim(0, 1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,family='serif'),
        axis.text.x = element_text(size=13,family='serif',
                                   angle = 45,hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=13,family='serif'),
        axis.title.y = element_text(size=13,family='serif'))



ggplot(data, mapping = aes(x = algorithm, y = rand, fill = methods)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#4682B4", "#F0E442", "#E69F00")) +
  ylab('RI') + xlab("") + theme_bw() + ylim(0, 1) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12,family='serif'),
        axis.text.x = element_text(size=13,family='serif',
                                   angle = 45,hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=13,family='serif'),
        axis.title.y = element_text(size=13,family='serif'))
