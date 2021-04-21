library(ggplot2)
library(igraph)
library(MASS)
library(pROC)
library(ROCR)
library(ggplot2)
source('help_func.R')
library(Matrix)




#----------------------------------------------------------------------
# Compare the result of NR-F and NE when N=250
#----------------------------------------------------------------------
n1 <- 50
l <- 5
loops <- 50
n <- n1*l
p_in <- seq(0.25,0.49,0.01)
ref_mode <- 1-diag(1,nrow = n1)
ref <- as.matrix(bdiag(ref_mode,ref_mode,ref_mode,ref_mode,ref_mode))

AUROC_NR <- matrix(0,nrow=length(p_in),ncol=loops)
AUPR_NR <- matrix(0,nrow=length(p_in),ncol=loops)
AUROC_NE <- matrix(0,nrow=length(p_in),ncol=loops)
AUPR_NE <- matrix(0,nrow=length(p_in),ncol=loops)

for (i in 1:length(p_in)){
  for (j in 1:loops){
    print(paste0('p_in: ',i,'th','   loops: ',j))
    
    A <- read.csv(paste0('./simudata_250/noisy_',i,'_loop_',j,'.csv'))
    A <- as.matrix(A[,-1])
    
    
    W_NR <- NR_F(A,m=2)
    W_NR <- (W_NR+t(W_NR))/2
    W_NR <- W_NR-diag(diag(W_NR))
    W_NR <- W_NR/max(W_NR)
  
    W_NE <- as.matrix(read.table(
      paste0('./simudata_250_after/NE',i,'_loops',j,'.csv'),sep = ','))
    W_NE <- (W_NE+t(W_NE))/2
    W_NE <- W_NE-diag(diag(W_NE))
    W_NE <- W_NE/max(W_NE)
    

    test <- data.frame(label=as.factor(ref*A),
                       output_NR=as.vector(W_NR),
                       output_ND=as.vector(W_NE))
    
    pred_NR <-prediction(test$output_NR,test$label)
    pred_NE <-prediction(test$output_ND,test$label)
    
    auc_NR <- performance(pred_NR,'auc')
    AUROC_NR[i,j] <- unlist(slot(auc_NR,"y.values"))
    auc_NE <- performance(pred_NE 
                          ,'auc')
    AUROC_NE[i,j] <- unlist(slot(auc_NE,"y.values"))
    
    pr_NR <- performance(pred_NR,'aucpr')
    AUPR_NR[i,j] <- unlist(slot(pr_NR,"y.values"))
    pr_NE <- performance(pred_NE,'aucpr')
    AUPR_NE[i,j] <- unlist(slot(pr_NE,"y.values"))

  }
}






#------------------------------------------------------
data <-data.frame(p_in=c(rep(p_in,2)),
                  group=c(rep('NE',25),rep('NR-F',25)),
                  score_auc=c(apply(AUROC_NE,1,mean),
                              apply(AUROC_NR,1,mean)),
                  sd_auc=c(apply(AUROC_NE,1,sd),
                           apply(AUROC_NR,1,sd)),
                  score_pr=c(apply(AUPR_NE,1,mean),
                             apply(AUPR_NR,1,mean)),
                  sd_pr=c(apply(AUPR_NE,1,sd),
                          apply(AUPR_NR,1,sd)))
#-----------------------------------------------------------------
# plot performance
ggplot(data,aes(x=p_in,y=score_auc,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_auc-sd_auc,ymax=score_auc+sd_auc), alpha = 0.05)+
  scale_fill_manual(values = c('#56B4E9', '#F0E442'))+
  scale_color_manual(values = c('#56B4E9', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUROC (l = 5, n = 50, N = 250)'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))


ggplot(data,aes(x=p_in,y=score_pr,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_pr-sd_pr,ymax=score_pr+sd_pr), alpha = 0.05)+
  scale_fill_manual(values = c('#56B4E9', '#F0E442'))+
  scale_color_manual(values = c('#56B4E9', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUPR (l = 5, n = 50, N = 250)'))+
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
library(Matrix)

#----------------------------------------------------------------------
# Compare the result of NR-F and NE when N=500
#----------------------------------------------------------------------
n1 <- 100
l <- 5
loops <- 50
n <- n1*l
p_in <- seq(0.25,0.49,0.01)
ref_mode <- 1-diag(1,nrow = n1)
ref <- as.matrix(bdiag(ref_mode,ref_mode,ref_mode,ref_mode,ref_mode))

AUROC_NR <- matrix(0,nrow=length(p_in),ncol=loops)
AUPR_NR <- matrix(0,nrow=length(p_in),ncol=loops)
AUROC_NE <- matrix(0,nrow=length(p_in),ncol=loops)
AUPR_NE <- matrix(0,nrow=length(p_in),ncol=loops)

for (i in 1:length(p_in)){
  for (j in 1:loops){
    print(paste0('p_in: ',i,'th','   loops: ',j))
    
    A <- read.csv(paste0('./simudata_500/noisy_',i,'_loop_',j,'.csv'))
    A <- as.matrix(A[,-1])
    
    
    W_NR <- NR_F(A,m=2)
    W_NR <- (W_NR+t(W_NR))/2
    W_NR <- W_NR-diag(diag(W_NR))
    W_NR <- W_NR/max(W_NR)
    
    W_NE <- as.matrix(read.table(
      paste0('./simudata_500_after/NE',i,'_loops',j,'.csv'),sep = ','))
    W_NE <- (W_NE+t(W_NE))/2
    W_NE <- W_NE-diag(diag(W_NE))
    W_NE <- W_NE/max(W_NE)
    
    
    test <- data.frame(label=as.factor(ref*A),
                       output_NR=as.vector(W_NR),
                       output_ND=as.vector(W_NE))
    
    pred_NR <-prediction(test$output_NR,test$label)
    pred_NE <-prediction(test$output_ND,test$label)
    
    auc_NR <- performance(pred_NR,'auc')
    AUROC_NR[i,j] <- unlist(slot(auc_NR,"y.values"))
    auc_NE <- performance(pred_NE 
                          ,'auc')
    AUROC_NE[i,j] <- unlist(slot(auc_NE,"y.values"))
    
    pr_NR <- performance(pred_NR,'aucpr')
    AUPR_NR[i,j] <- unlist(slot(pr_NR,"y.values"))
    pr_NE <- performance(pred_NE,'aucpr')
    AUPR_NE[i,j] <- unlist(slot(pr_NE,"y.values"))
    
  }
}



#------------------------------------------------------
data <-data.frame(p_in=c(rep(p_in,2)),
                  group=c(rep('NE',25),rep('NR-F',25)),
                  score_auc=c(apply(AUROC_NE,1,mean),
                              apply(AUROC_NR,1,mean)),
                  sd_auc=c(apply(AUROC_NE,1,sd),
                           apply(AUROC_NR,1,sd)),
                  score_pr=c(apply(AUPR_NE,1,mean),
                             apply(AUPR_NR,1,mean)),
                  sd_pr=c(apply(AUPR_NE,1,sd),
                          apply(AUPR_NR,1,sd)))

#-----------------------------------------------------------------
# plot performance
ggplot(data,aes(x=p_in,y=score_auc,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_auc-sd_auc,ymax=score_auc+sd_auc), alpha = 0.05)+
  scale_fill_manual(values = c('#56B4E9', '#F0E442'))+
  scale_color_manual(values = c('#56B4E9', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUROC (l = 5, n = 100, N = 500)'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))


ggplot(data,aes(x=p_in,y=score_pr,
                group=group,color=group))+
  geom_line(size=0.8)+geom_point(size=1.5)+
  geom_ribbon(aes(ymin=score_pr-sd_pr,ymax=score_pr+sd_pr), alpha = 0.05)+
  scale_fill_manual(values = c('#56B4E9', '#F0E442'))+
  scale_color_manual(values = c('#56B4E9', '#F0E442'))+
  theme_bw()+
  labs(title = paste0('AUPR (l = 5, n = 100, N = 500)'))+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=15,family='serif'))












