library(ggplot2)
library(patchwork)
library(RColorBrewer)

ssim <-t(read.csv(
      "./ssim.csv",
      encoding = "UTF-8",
      header = F))
spear <-
  t(
    read.csv(
      "./spearman.csv",
      encoding = "UTF-8",
      header = F
    )
  )
pear <-
  t(
    read.csv(
      "./pearson.csv",
      encoding = "UTF-8",
      header = F
    )
  )


data <- data.frame(
  methods = c(
    rep('downsampled', 22),
    rep('NR-F', 22),
    rep('NR-F*2', 22),
    rep('NE', 22),
    rep('ND', 22)
  ),
  SSIM = as.vector(t(ssim)),
  spearman = as.vector(t(spear)),
  pearson = as.vector(t(pear))
)

data$methods <- factor(data$methods,
                       c('downsampled', 'ND',
                         'NE',
                         'NR-F',
                         'NR-F*2'))

p1 <- ggplot(data, mapping = aes(methods,SSIM)) +
  stat_boxplot(geom = "errorbar", width = 0.8) +
  geom_boxplot(aes(fill = methods),outlier.colour = NA,
    size = 1,width=1,
    fill = "white"
  ) +
  geom_jitter(aes(fill = methods,color=methods),
              width = 0.2,shape = 21,
              size = 1.5)+ 
  geom_hline(aes(yintercept=median(ssim[1,])),colour="#990000", linetype="dashed",size=1)+
  scale_fill_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
  scale_color_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position='none',
        title =element_text(size=8, family='serif'),
        axis.text.x = element_text(
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,size=10,family='serif'),
    axis.text.y = element_text(size=10,family='serif')) + 
  guides(fill = F)+xlab('')+ylab('')+labs(title = "SSIM")+ylim(0,1)


p2 <- ggplot(data, mapping = aes(methods,spearman)) +
  stat_boxplot(geom = "errorbar", width = 0.8) +
  geom_boxplot(aes(fill = methods),outlier.colour = NA,
               size = 0.8,width=1,
               fill = "white"
  ) +
  geom_jitter(aes(fill = methods,color=methods),
              width = 0.2,shape = 21,
              size = 1.5)+ 
  geom_hline(aes(yintercept=median(spear[1,])),colour="#990000", linetype="dashed",size=1)+
  scale_fill_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
  scale_color_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=8, family='serif'),
        legend.position='none',
        axis.text.x = element_text(
          angle = 90,
          hjust = 0.5,
          vjust = 0.5,size=10,family='serif'),
        axis.text.y = element_text(size=10,family='serif')) + 
  guides(fill = F)+xlab('')+ylab('')+labs(title = "Spearman correlation")+ylim(0,1)
  
  
  p3 <- ggplot(data, mapping = aes(methods,pearson)) +
    stat_boxplot(geom = "errorbar", width = 0.8) +
    geom_boxplot(aes(fill = methods),outlier.colour = NA,
                 size = 0.8,width=1,
                 fill = "white"
    ) +
    geom_jitter(aes(fill = methods,color=methods),
                width = 0.2,shape = 21,
                size = 1.5)+ 
    geom_hline(aes(yintercept=median(pear[1,])),colour="#990000", linetype="dashed",size=1)+
    scale_fill_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
    scale_color_manual(values = c("#999999", "#009E73", "#56B4E9","#F0E442", "#E69F00"))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position='none',
          title =element_text(size=8, family='serif'),
          axis.text.x = element_text(
            angle = 90,
            hjust = 0.5,
            vjust = 0.5,size=10,family='serif'),
          axis.text.y = element_text(size=10,family='serif')) + 
    guides(fill = F)+xlab('')+ylab('')+
    labs(title = "Pearson correlation")+ylim(0,1)
  
  
  p1 + p2 + p3
