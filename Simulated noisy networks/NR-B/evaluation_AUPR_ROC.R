library(ggplot2)
library(igraph)
library(MASS)
library(pROC)
library(ROCR)
library(ggplot2)
source('help_func.R')

#----------------------------------------------------------------------
# 1.Compare the result of NR-B and ND on noisy BA graph
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Average the result of 100 BA graph
n <- 50
loops <- 100

noisy_proportion <- seq(0.1,0.5,0.01)#noise proportion
auc_NR <- matrix(0,nrow=loops,
                    ncol=length(noisy_proportion))
pr_NR <- matrix(0,nrow=loops,
                   ncol=length(noisy_proportion))
auc_ND <- matrix(0,nrow=loops,
                    ncol=length(noisy_proportion))
pr_ND <- matrix(0,nrow=loops,
                   ncol=length(noisy_proportion))


# import BA data
for (i in 1:length(noisy_proportion)){
  for (j in 1:loops){
    print(paste0('i_',i,'    j_',j))
    A1 <- read.csv(paste0('./simudata_BA/true_',i,'_loop_',j,'.csv'))
    A2 <- read.csv(paste0('./simudata_BA/noisy_',i,'_loop_',j,'.csv'))
    A1 <- as.matrix(A1[,2:51])
    A2 <- as.matrix(A2[,2:51])

    #NR-B
    print('****************************************************')
    pre <- pre_processing(A2)
    P <- pre[[1]]
    print(c(pre[[2]],pre[[3]]))
    print('****************************************************')
    D <- as.vector(Null(P-diag(1,nrow=n)))
    D <- D/sum(D)
    W_NR <- diag(D)%*%P
    W_NR <- (W_NR+t(W_NR))/2
    W_NR <- W_NR-diag(diag(W_NR))
    W_NR <- W_NR/max(W_NR)

    #ND
    W_ND <- as.matrix(read.table(paste0('./simudata_BA_after/ND_noisy',
            i,'_loop',j,'.txt'),sep=','))

    test <- data.frame(label=as.factor(as.matrix(A1)),
                       output_NR=as.vector(W_NR),
                       output_ND=as.vector(W_ND))

    pred_NR <-prediction(test$output_NR,test$label)
    pred_ND <-prediction(test$output_ND,test$label)
    auc_NR1 <- performance(pred_NR,'auc')
    auc_NR[j,i] <- unlist(slot(auc_NR1,"y.values"))
    auc_ND1 <- performance(pred_ND,'auc')
    auc_ND[j,i] <- unlist(slot(auc_ND1,"y.values"))
    pr_NR1 <- performance(pred_NR,'aucpr')
    pr_NR[j,i] <- unlist(slot(pr_NR1,"y.values"))
    pr_ND1 <- performance(pred_ND,'aucpr')
    pr_ND[j,i] <- unlist(slot(pr_ND1,"y.values"))


  }
}


#------------------------------------------------------
data <-data.frame(noise_proportion=c(rep(noisy_proportion,2)),
                  group=c(rep('ND',41),rep('NR-B',41)),
                  score_auc=c(apply(auc_ND,2,mean),
                          apply(auc_NR,2,mean)),
                  sd_auc=c(apply(auc_ND,2,sd),
                       apply(auc_NR,2,sd)),
                  score_pr=c(apply(pr_ND,2,mean),
                              apply(pr_NR,2,mean)),
                  sd_pr=c(apply(pr_ND,2,sd),
                           apply(pr_NR,2,sd)))
#-----------------------------------------------------------------
# plot performance
ggplot(data,aes(x=noise_proportion,y=score_auc,
                   group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_auc-sd_auc,ymax=score_auc+sd_auc), alpha = 0.05)+
  scale_fill_manual(values = c('#009E73', '#F0E442'))+
  scale_color_manual(values = c('#009E73', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUROC of BA graphs'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))


ggplot(data,aes(x=noise_proportion,y=score_pr,
                   group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_pr-sd_pr,ymax=score_pr+sd_pr), alpha = 0.05)+
  scale_fill_manual(values = c('#009E73', '#F0E442'))+
  scale_color_manual(values = c('#009E73', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUPR of BA graphs'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
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

auc_NR <- matrix(0,nrow=loops,
                 ncol=length(noisy_proportion))
pr_NR <- matrix(0,nrow=loops,
                ncol=length(noisy_proportion))
auc_ND <- matrix(0,nrow=loops,
                 ncol=length(noisy_proportion))
pr_ND <- matrix(0,nrow=loops,
                ncol=length(noisy_proportion))

# import ER graphs
for (i in 1:length(noisy_proportion)){
  for (j in 1:loops){
    print(paste0('i_',i,'    j_',j))
    A1 <- read.csv(paste0('./simudata_ER/true_',i,'_loop_',j,'.csv'))
    A2 <- read.csv(paste0('./simudata_ER/noisy_',i,'_loop_',j,'.csv'))
    A1 <- as.matrix(A1[,2:51])
    A2 <- as.matrix(A2[,2:51])
    
    #NR-B
    print('****************************************************')
    pre <- pre_processing(A2)
    P <- pre[[1]]
    print('****************************************************')
    D <- as.vector(Null(P-diag(1,nrow=n)))
    D <- D/sum(D)
    W_NR <- diag(D)%*%P
    W_NR <- (W_NR+t(W_NR))/2
    W_NR <- W_NR-diag(diag(W_NR))
    W_NR <- W_NR/max(W_NR)
    
    #ND
    W_ND <- as.matrix(read.table(paste0('./simudata_ER_after/ND_noisy',
                                        i,'_loop',j,'.txt'),sep=','))
    
    test <- data.frame(label=as.factor(as.matrix(A1)),
                       output_NR=as.vector(W_NR),
                       output_ND=as.vector(W_ND))
    
    pred_NR <-prediction(test$output_NR,test$label)
    pred_ND <-prediction(test$output_ND,test$label)
    auc_NR1 <- performance(pred_NR,'auc')
    auc_NR[j,i] <- unlist(slot(auc_NR1,"y.values"))
    auc_ND1 <- performance(pred_ND,'auc')
    auc_ND[j,i] <- unlist(slot(auc_ND1,"y.values"))
    pr_NR1 <- performance(pred_NR,'aucpr')
    pr_NR[j,i] <- unlist(slot(pr_NR1,"y.values"))
    pr_ND1 <- performance(pred_ND,'aucpr')
    pr_ND[j,i] <- unlist(slot(pr_ND1,"y.values"))
    
    
  }
}


#------------------------------------------------------
data_ER <-data.frame(noise_proportion=c(rep(noisy_proportion,2)),
                     group=c(rep('ND',41),rep('NR-B',41)),
                     score_auc=c(apply(auc_ND,2,mean),
                                 apply(auc_NR,2,mean)),
                     sd_auc=c(apply(auc_ND,2,sd),
                              apply(auc_NR,2,sd)),
                     score_pr=c(apply(pr_ND,2,mean),
                                apply(pr_NR,2,mean)),
                     sd_pr=c(apply(pr_ND,2,sd),
                             apply(pr_NR,2,sd)))
#-----------------------------------------------------------------
# plot performance
ggplot(data_ER,aes(x=noise_proportion,y=score_auc,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_auc-sd_auc,ymax=score_auc+sd_auc), alpha = 0.05)+
  scale_fill_manual(values = c('#009E73', '#F0E442'))+
  scale_color_manual(values = c('#009E73', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUROC of ER graphs'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))


ggplot(data_ER,aes(x=noise_proportion,y=score_pr,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_pr-sd_pr,ymax=score_pr+sd_pr), alpha = 0.05)+
  scale_fill_manual(values = c('#009E73', '#F0E442'))+
  scale_color_manual(values = c('#009E73', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUPR of ER graphs'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))
