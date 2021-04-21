source('help_func.R')


#------------------------------------------------------------------
# generate noisy graphs (N=250)
#------------------------------------------------------------------


n1 <- 50
loops <- 50
l <- 5
n <- n1*l
p_in <- seq(0.25,0.49,0.01)
p_out <- (0.5-p_in)/(l-1)

for (i in 1:length(p_in)){
  for (j in 1:loops){
    print(paste0('l=',l,'   p_in: ',i,'th','   loops: ',j))
    A <- generate_plp(n1,l,p_in[i],p_out[i])
    write.csv(A,paste0('./simudata_250/noisy_',i,'_loop_',j,'.csv'))

  }
}
