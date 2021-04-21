library(ROCR)
library(MASS)
library(ggplot2)
source('NR_B.R')

#----------------------------------------------------------------------------
# Application3: 
# NR-B distinguishes strong collaborations of co-authorship networks
# We calculate the AURUC and AUPR to evaluate the performance of NR-B and NE
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# import weighted co-authorship matrix and 
# create input unweighted co-authorship matrix
net_ref_3col <- read.csv('coauthor.csv', header = F)
n <- max(net_ref_3col) + 1

# weighted
net_ref_weighted <- matrix(0, nrow = n, ncol = n)
for (i in 1:dim(net_ref_3col)[1]) {
  net_ref_weighted[net_ref_3col[i, 1] + 1, net_ref_3col[i, 2] + 1] =
    net_ref_3col[i, 3]
}
net_ref_weighted = net_ref_weighted - diag(diag(net_ref_weighted))
net_ref_weighted = net_ref_weighted + t(net_ref_weighted)
rs <- rowSums(net_ref_weighted)
cs <- colSums(net_ref_weighted)
r_index <- which(rs == 0)
c_index <- which(cs == 0)

net_ref_weighted <- net_ref_weighted[-r_index, ]
net_ref_weighted <- net_ref_weighted[, -c_index]

# unweighted
net_ref_unweighted <- (net_ref_weighted != 0) * 1


#-----------------------------------------------------------------------------
# Applying NR-B algorithm and importing ND result
net_denoise_NR_B <- NR_B(net_ref_unweighted, 2, 0.1, 1)
net_denoise_NR_B <- net_denoise_NR_B * net_ref_unweighted

net_denoise_ND <- as.matrix(read.csv('ND.csv', header = F))


#-------------------------------------------------------
# compute the strong/weak relationships
net_ref_strong <- (net_ref_weighted >= 1 / 2) * 1
net_ref_weak <- (0 < net_ref_weighted & net_ref_weighted < 1 / 2) * 1

n_strong <- sum(net_ref_strong) / 2
n_weak <- sum(net_ref_weak) / 2
n_edges <- sum(net_ref_unweighted) / 2

eva_ref_net <- (net_ref_strong == 1) * 2 + (net_ref_weak == 1) * 1
eva_ref_vec <- eva_ref_net[upper.tri(eva_ref_net, diag = F)]
eva_ref_vec <- eva_ref_vec[eva_ref_vec > 0]

# normalization for NR-B
eva_denoise_NR_B_vec <-
  net_denoise_NR_B[upper.tri(net_denoise_NR_B, diag = F)]
eva_denoise_NR_B_vec <-
  eva_denoise_NR_B_vec[eva_denoise_NR_B_vec > 0]
eva_denoise_NR_B_vec <-
  (eva_denoise_NR_B_vec - min(eva_denoise_NR_B_vec)) /
  (max(eva_denoise_NR_B_vec) - min(eva_denoise_NR_B_vec))

# normalization for ND
eva_denoise_ND_vec <-
  net_denoise_ND[upper.tri(net_denoise_ND, diag = F)]
eva_denoise_ND_vec <-
  eva_denoise_ND_vec[eva_denoise_ND_vec > 0]
eva_denoise_ND_vec <-
  (eva_denoise_ND_vec - min(eva_denoise_ND_vec)) /
  (max(eva_denoise_ND_vec) - min(eva_denoise_ND_vec))


eva_data <- data.frame(ref = eva_ref_vec - 1,
                       denoise_NR_B = eva_denoise_NR_B_vec,
                       denoise_ND = eva_denoise_ND_vec)


#--------------------------------------------------------------------
# evaluation of NR-B
print('Evaluation of AUROC and AUPR for NR-B.....................')

pred_NR_B <- prediction(eva_data$denoise_NR_B, eva_data$ref)

perf_auc_NR_B <- performance(pred_NR_B, "tpr", "fpr")
auc_NR_B <- performance(pred_NR_B, 'auc')
auc_NR_B <- unlist(slot(auc_NR_B, "y.values"))
print(paste0('AUROC of NR-B:  ', auc_NR_B))

perf_pr_NR_B <- performance(pred_NR_B, "prec", "tpr")
pr_NR_B <- performance(pred_NR_B, 'aucpr')
pr_NR_B <- unlist(slot(pr_NR_B, "y.values"))
print(paste0('AUPR of NR-B:  ', pr_NR_B))



#--------------------------------------------------------
# evaluation of ND
print('Evaluation of AUROC and AUPR for ND.....................')

pred_ND <- prediction(eva_data$denoise_ND, eva_data$ref)

perf_auc_ND <- performance(pred_ND, "tpr", "fpr")
auc_ND <- performance(pred_ND, 'auc')
auc_ND <- unlist(slot(auc_ND, "y.values"))
print(paste0('AUROC of ND:  ', auc_ND))

perf_pr_ND <- performance(pred_ND, "prec", "tpr")
pr_ND <- performance(pred_ND, 'aucpr')
pr_ND <- unlist(slot(pr_ND, "y.values"))
print(paste0('AUPR of ND:  ', pr_ND))





#-------------------------------------------------------------
# plot AUC curves
data_auc <-
  data.frame(
    methods = as.factor(c(rep(
      'ND', length(perf_auc_ND@x.values[[1]])
    ),
    rep(
      'NR-B', length(perf_auc_NR_B@x.values[[1]])
    ))),
    x = c(perf_auc_ND@x.values[[1]],
          perf_auc_NR_B@x.values[[1]]),
    y = c(perf_auc_ND@y.values[[1]],
          perf_auc_NR_B@y.values[[1]])
  )

ggplot(data_auc, aes(x, y, group = methods, color = methods)) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  geom_line(size = 1.2) +
  geom_segment(aes(
    x = 0,
    xend = 1,
    y = 0,
    yend = 1
  ),
  colour = 'black',
  linetype = "dashed") +
  labs(title = "AUC curve", x = "FPR", y = "TPR") +
  annotate(
    "text",
    x = 0.2,
    y = 0.6,
    label = paste0("AUROC:", round(auc_NR_B, 5)),
    parse = T,
    color = '#E69F00'
  ) +
  annotate(
    "text",
    x = 0.25,
    y = 0.9,
    label = paste0("AUROC:", round(auc_ND, 5)),
    parse = T,
    color = '#009E73'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 15,
      family = 'serif',
      hjust = 0.5
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, family = 'serif'),
    axis.text.x = element_text(size = 15, family = 'serif'),
    axis.text.y = element_text(size = 15, family = 'serif'),
    axis.title.y = element_text(size = 15, family = 'serif'),
    axis.title.x = element_text(size = 15, family = 'serif')
  )



#-------------------------------------------------------------------
# plot PR curves
data_pr <-
  data.frame(
    methods = as.factor(c(rep(
      'ND', length(perf_pr_ND@x.values[[1]])
    ),
    rep(
      'NR-B', length(perf_pr_NR_B@x.values[[1]])
    ))),
    x = c(perf_pr_ND@x.values[[1]],
          perf_pr_NR_B@x.values[[1]]),
    y = c(perf_pr_ND@y.values[[1]],
          perf_pr_NR_B@y.values[[1]])
  )
data_pr <- na.omit(data_pr)

ggplot(data_pr, aes(x, y, group = methods, color = methods)) +
  scale_fill_manual(values = c("#009E73", "#E69F00")) +
  scale_color_manual(values = c("#009E73", "#E69F00")) +
  geom_line(size = 1.2) +
  labs(title = "PR curve", x = "TPR", y = "Precision") +
  annotate(
    "text",
    x = 0.9,
    y = 1.0,
    label = paste0("AUPR:", round(pr_NR_B, 5)),
    parse = T,
    color = '#E69F00'
  ) +
  annotate(
    "text",
    x = 0.5,
    y = 0.7,
    label = paste0("AUPR:", round(pr_ND, 5)),
    parse = T,
    color = '#009E73'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      size = 15,
      family = 'serif',
      hjust = 0.5
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 15, family = 'serif'),
    axis.text.x = element_text(size = 15, family = 'serif'),
    axis.text.y = element_text(size = 15, family = 'serif'),
    axis.title.y = element_text(size = 15, family = 'serif'),
    axis.title.x = element_text(size = 15, family = 'serif')
  )
