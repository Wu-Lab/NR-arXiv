library(ggplot2)
library(igraph)
library(MASS)
source('help_func.R')


#----------------------------------------------------------------------
# Because this runs slowly, 
# we strongly recommend that you reduce the number of loops or 
# just load the result we got and go to line 90 to plot the performance
load("./results/N=250.RData")



#----------------------------------------------------------------------
# Compare the result of NR-F and NE when N=250
#----------------------------------------------------------------------

n1 <- 50
l <- 5
loops <- 50
n <- n1*l
p_in <- seq(0.25,0.49,0.01)
p_out <- (0.5-p_in)/(l-1)


error_before <- matrix(0,nrow=length(p_in),ncol = loops)
error_after_NR2 <- matrix(0,nrow=length(p_in),ncol = loops)
error_after_NR3<- matrix(0,nrow=length(p_in),ncol = loops)
error_after_NE <- matrix(0,nrow=length(p_in),ncol = loops)

for (i in 1:length(p_in)){
  for (j in 1:loops){
    print(paste0('p_in: ',i,'th','   loops: ',j))
    A <- read.csv(paste0('./simudata_250/noisy_',i,'_loop_',j,'.csv'))
    A <- as.matrix(A[,-1])
    
    error_before[i,j] <- compute_error(n1,l,A)
    
    
    output_w_NR2 <- NR_F(A,2)
    output_w_NR2 <- (output_w_NR2+t(output_w_NR2))/2
    output_w_NR2 <- output_w_NR2-diag(diag(output_w_NR2))
    output_w_NR2 <- output_w_NR2/max(output_w_NR2)
    output_adj_NR2 <- find_thre(output_w_NR2)
    
    output_w_NR3 <- NR_F(A,3)
    output_w_NR3 <- (output_w_NR3+t(output_w_NR3))/2
    output_w_NR3 <- output_w_NR3-diag(diag(output_w_NR3))
    output_w_NR3 <- output_w_NR3/max(output_w_NR3)
    output_adj_NR3 <- find_thre(output_w_NR3)
    
    
    
    output_w_NE <- read.table(paste0('./simudata_250_after/NE',i,'_loops',j,'.csv'),sep = ',')
    output_w_NE <- as.matrix(output_w_NE)
    output_w_NE <- (output_w_NE+t(output_w_NE))/2
    output_w_NE <- output_w_NE-diag(diag(output_w_NE))
    output_w_NE <- output_w_NE/max(output_w_NE)
    output_adj_NE <- find_thre(output_w_NE)
    
    error_after_NR2[i,j] <- compute_error(n1,l,output_adj_NR2)
    error_after_NR3[i,j] <- compute_error(n1,l,output_adj_NR3)
    error_after_NE[i,j] <- compute_error(n1,l,output_adj_NE)
    
  }
}


data <-data.frame(p_in=c(rep(p_in,4)),
                  group=c(rep('before',25),
                          rep('NR-F(m=2)',25),
                          rep('NR-F(m=3)',25),
                          rep('NE',25)),
                  score=c(apply(error_before,1,mean),
                          apply(error_after_NR2,1,mean),
                          apply(error_after_NR3,1,mean),
                          apply(error_after_NE,1,mean)),
                  sd=c(apply(error_before,1,sd)/2,
                       apply(error_after_NR2,1,sd)/2,
                       apply(error_after_NR3,1,sd)/2,
                       apply(error_after_NE,1,sd)/2))





#---------------------------------------------------------------------
# plot performanc
ggplot(data,aes(x=p_in,y=score,
                group=group,color=group))+
  geom_line(size=1)+geom_point(size=2)+
  geom_ribbon(aes(ymin=score-sd/2,
                  ymax=score+sd/2), alpha = 0.05)+
  scale_fill_manual(values = c("#999999", 
                               "#56B4E9",
                               "#F0E442", 
                               "#E69F00"))+
  scale_color_manual(values = c("#999999", 
                                "#56B4E9",
                                "#F0E442", 
                                "#E69F00"))+
  theme_bw()+
  labs(title = paste0('l = ',l,', n = ',n1,', N = ',n))+
  ylab('1/SNR')+
  theme(plot.title = element_text(size=15,family='serif',hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=14,family='serif'),
        axis.text.x = element_text(size=15,family='serif'),
        axis.text.y = element_text(size=15,family='serif'),
        axis.title.y = element_text(size=15,family='serif'),
        axis.title.x = element_text(size=15,family='serif'))


