source('help_func.R')
library(ggplot2)
library(igraph)
library(MASS)
library(pROC)
library(ROCR)
library(ggplot2)

#----------------------------------------------------------------------
# 1.Compare the result of NR-B and ND on noisy BA graph
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Average the result of 100 BA graph
n <- 50
loops <- 100

noisy_proportion <- seq(0.1,0.5,0.01)#noise proportion 
Error_before <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_NR_B2 <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_NR_B3 <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_ND <- matrix(0,nrow=length(noisy_proportion),ncol=loops)


# import BA data
for (i in 1:length(noisy_proportion)){
  for (j in 1:loops){
    print(paste0('i_',i,'    j_',j))
    A1 <- read.csv(paste0('./simudata_BA/true_',i,'_loop_',j,'.csv'))
    A2 <- read.csv(paste0('./simudata_BA/noisy_',i,'_loop_',j,'.csv'))
    A1 <- as.matrix(A1[,2:51])
    A2 <- as.matrix(A2[,2:51])
    
    edges = sum(A1)/2
    
    #NR-B(m=2)
    print('****************************************************')
    pre2 <- pre_processing(A2)
    P2 <- pre2[[1]]
    print(c(pre2[[2]],pre2[[3]]))
    print('****************************************************')
    D2 <- as.vector(Null(P2-diag(1,nrow=n)))
    D2 <- D2/sum(D2)
    W_NR2 <- diag(D2)%*%P2
    W_NR2 <- (W_NR2+t(W_NR2))/2
    W_NR2 <- W_NR2-diag(diag(W_NR2))
    W_NR2 <- W_NR2/max(W_NR2)
    
    #NR-B(m=3)
    print('****************************************************')
    pre3 <- pre_processing(A2,m=3)
    P3 <- pre3[[1]]
    print(c(pre3[[2]],pre3[[3]]))
    print('****************************************************')
    D3 <- as.vector(Null(P3-diag(1,nrow=n)))
    D3 <- D3/sum(D3)
    W_NR3 <- diag(D3)%*%P3
    W_NR3 <- (W_NR3+t(W_NR3))/2
    W_NR3 <- W_NR3-diag(diag(W_NR3))
    W_NR3 <- W_NR3/max(W_NR3)
    
    
    #ND
    W_ND <- as.matrix(read.table(paste0('./simudata_BA_after/ND_noisy',
                                        i,'_loop',j,'.txt'),sep=','))
    
    adj_NR_B2 <- change_into_adj_keep_edges(W_NR2,edges)
    adj_NR_B3 <- change_into_adj_keep_edges(W_NR3,edges)
    adj_ND <- change_into_adj_keep_edges(W_ND,edges)
    
    Error_before[i,j] = compute_error(A1,A2)
    Error_NR_B2[i,j] = compute_error(A1,adj_NR_B2)
    Error_NR_B3[i,j] = compute_error(A1,adj_NR_B3)
    Error_ND[i,j] = compute_error(A1,adj_ND)

    
  }
}


#-----------------------------------------------------------------
data <-data.frame(noise=c(rep(noisy_proportion,4)),
                  group=c(rep('before',41),
                          rep('NR-B(m=2)',41),
                          rep('NR-B(m=3)',41),
                          rep('ND',41)),
                  score=c(apply(Error_before,1,mean),
                          apply(Error_NR_B2,1,mean),
                          apply(Error_NR_B3,1,mean),
                          apply(Error_ND,1,mean)),
                  sd=c(apply(Error_before,1,sd)/2,
                       apply(Error_NR_B2,1,sd)/2,
                       apply(Error_NR_B3,1,sd)/2,
                       apply(Error_ND,1,sd)/2))

data$group <- factor(data$group,levels = c('before',
                                           'NR-B(m=2)',
                                           'NR-B(m=3)',
                                           'ND'))

ggplot(data,aes(x=noise,y=score,group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin = score-sd, ymax=score+sd), alpha = 0.05)+
  scale_fill_manual(values = c("#999999",
                               '#F0E442',
                               '#E69F00',
                               '#009E73'))+
  scale_color_manual(values = c("#999999",
                                '#F0E442',
                                '#E69F00',
                                '#009E73'))+
  theme_bw()+
  labs(title = paste0('BA graphs (N = 50, m = 3)'))+
  xlab('noise_proportion')+
  ylab('1/SNR')+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_text(size=15,family='serif'),
        axis.title.x = element_text(size=15,family='serif'))












library(ggplot2)
library(igraph)
library(MASS)
library(pROC)
library(ROCR)
library(ggplot2)
source('help_func.R')
#----------------------------------------------------------------------
#     2.Compare the result of NR-B and ND on noisy ER graph
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Average the result of 100 ER random graphs
p <- 0.3
n <- 50
loops <- 100
noisy_proportion <- seq(0.1,0.5,0.01)#noise proportion 
Error_before <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_NR_B2 <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_NR_B3 <- matrix(0,nrow=length(noisy_proportion),ncol=loops)
Error_ND <- matrix(0,nrow=length(noisy_proportion),ncol=loops)

# import ER graphs
for (i in 1:length(noisy_proportion)){
  for (j in 1:loops){
    print(paste0('i_',i,'    j_',j))
    A1 <- read.csv(paste0('./simudata_ER/true_',i,'_loop_',j,'.csv'))
    A2 <- read.csv(paste0('./simudata_ER/noisy_',i,'_loop_',j,'.csv'))
    A1 <- as.matrix(A1[,2:51])
    A2 <- as.matrix(A2[,2:51])
    
    edges = sum(A1)/2
    
    #NR-B(m=2)
    print('****************************************************')
    pre2 <- pre_processing(A2)
    P2 <- pre2[[1]]
    print(c(pre2[[2]],pre2[[3]]))
    print('****************************************************')
    D2 <- as.vector(Null(P2-diag(1,nrow=n)))
    D2 <- D2/sum(D2)
    W_NR2 <- diag(D2)%*%P2
    W_NR2 <- (W_NR2+t(W_NR2))/2
    W_NR2 <- W_NR2-diag(diag(W_NR2))
    W_NR2 <- W_NR2/max(W_NR2)
    
    #NR-B(m=3)
    print('****************************************************')
    pre3 <- pre_processing(A2,m=3)
    P3 <- pre3[[1]]
    print(c(pre3[[2]],pre3[[3]]))
    print('****************************************************')
    D3 <- as.vector(Null(P3-diag(1,nrow=n)))
    D3 <- D3/sum(D3)
    W_NR3 <- diag(D3)%*%P3
    W_NR3 <- (W_NR3+t(W_NR3))/2
    W_NR3 <- W_NR3-diag(diag(W_NR3))
    W_NR3 <- W_NR3/max(W_NR3)
    
    
    #ND
    W_ND <- as.matrix(read.table(paste0('./simudata_ER_after/ND_noisy',
                                        i,'_loop',j,'.txt'),sep=','))
    
    adj_NR_B2 <- change_into_adj_keep_edges(W_NR2,edges)
    adj_NR_B3 <- change_into_adj_keep_edges(W_NR3,edges)
    adj_ND <- change_into_adj_keep_edges(W_ND,edges)
    
    Error_before[i,j] = compute_error(A1,A2)
    Error_NR_B2[i,j] = compute_error(A1,adj_NR_B2)
    Error_NR_B3[i,j] = compute_error(A1,adj_NR_B3)
    Error_ND[i,j] = compute_error(A1,adj_ND)
    
    
    
  }
}

#-----------------------------------------------------------------
data_ER <-data.frame(noise=c(rep(noisy_proportion,4)),
                  group=c(rep('before',41),
                          rep('NR-B(m=2)',41),
                          rep('NR-B(m=3)',41),
                          rep('ND',41)),
                  score=c(apply(Error_before,1,mean),
                          apply(Error_NR_B2,1,mean),
                          apply(Error_NR_B3,1,mean),
                          apply(Error_ND,1,mean)),
                  sd=c(apply(Error_before,1,sd)/2,
                       apply(Error_NR_B2,1,sd)/2,
                       apply(Error_NR_B3,1,sd)/2,
                       apply(Error_ND,1,sd)/2))

data_ER$group <- factor(data_ER$group,levels = c('before',
                                           'NR-B(m=2)',
                                           'NR-B(m=3)',
                                           'ND'))

ggplot(data_ER,aes(x=noise,y=score,group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin = score-sd, ymax=score+sd), alpha = 0.05)+
  scale_fill_manual(values = c("#999999",
                               '#F0E442',
                               '#E69F00',
                               '#009E73'))+
  scale_color_manual(values = c("#999999",
                                '#F0E442',
                                '#E69F00',
                                '#009E73'))+
  theme_bw()+
  labs(title = paste0('ER graphs (N = 50, p = 0.3)'))+
  xlab('noise_proportion')+
  ylab('1/SNR')+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_text(size=15,family='serif'),
        axis.title.x = element_text(size=15,family='serif'))