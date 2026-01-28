rm(list=ls())

library(mvtnorm)
library(Matrix)
library(dplyr)

# work_path <- "/Users/xiaoyuexiong/XXYDATAanalysis/IP1-DetrendingPanelTimeSeries/CodesAndResultsForRevision/Everything2502/data"
# print(paste0("work_path: ", work_path))
# setwd(work_path)
work_path <- "D:/XXYDATAanalysis/IP_1/ClusterFile_PC/data"
print(paste0("work_path: ", work_path))
setwd(work_path)

######################## Control Parameters ########################
dateMark = "2502"
C = "dim1_mlGAR_v2"

nT <- nT_REPLACE  # 第一个参数是nT
nP <- nP_REPLACE  # 第二个参数是nP

ydim = 1 # number of dimensions of within-level variables
npad = 0  # initial number of time points to discard before retaining
N_repl= 100 # number of replications | replication k 
rs = 1:N_repl

######################## Model Parameters ########################
y = matrix(NA, nP, nT*2) #data matrix
#para_s = matrix(NA, nP, 4) #data matrix
e = matrix(NA, nP, nT*2) #residual data matrix
meer = matrix(NA, nP, nT) #data matrix

theta10 = 35 # asymptote, limitation when time gets close to +infinite
theta20 = 4 # sets the displacement along the x-axis (translates the graph to the left or right). 
theta30 = .8 #sets the growth rate (y scaling)
AR_mean = .3 # [mean of AR] 

r12 = 0 #cor(theta_1,theta_2)
r13 = 0 #cor(theta_1,theta_3)
r14 = 0.35 #cor(theta_1,AR)
r23 = 0 #cor(theta_2,theta_3)
r24 = 0 #cor(theta_2,AR)
r34 = 0 #cor(theta_3,AR)

sigma1 = 9 #SD(theta_1)
sigma2 = .5  #SD(theta_2)
sigma3 = .10  #SD(theta_3)
sigma4 = 0.1 #SD(AR)

Q = matrix(c(sigma1^2, r12*sigma1*sigma2, r13*sigma1*sigma3, r14*sigma1*sigma4,
             r12*sigma1*sigma2, sigma2^2, r23*sigma2*sigma3, r24*sigma2*sigma4,
             r13*sigma1*sigma3, r23*sigma2*sigma3, sigma3^2, r34*sigma3*sigma4,
             r14*sigma1*sigma4, r24*sigma2*sigma4, r34*sigma3*sigma4, sigma4^2), byrow=T,ncol=4)
Qr = chol(Q)

R = 1 #  process noise

########################################################################
######################## MC Experiments begins ########################
########################################################################
Total_t_begin <- proc.time()
for(k in 1:N_repl){
  
  r=rs[k]
  set.seed(r)
  
  dat_name <- paste0("Data_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark)
  
  dat_filename = paste0(dat_name,".csv")
  print(paste0("dat_filename: ", dat_filename))
  ######################## Data Generation ########################
  #------------------- Generate Data -------------------
  par(mfrow=c(1,1))
  for (p in 1:nP){
    time = seq(0.1,10,length=nT)
    up = rnorm(4)%*%Qr
    theta1 = theta10 + up[1]
    theta2 = theta20 + up[2]
    theta3 = theta30 + up[3]
    AR = AR_mean + up[4]
    #para_s[p,]=c(theta1,theta2,theta3,AR)
    for (i in 1:nT){
      tnow = time[i]
      mu_t = theta1*exp(-theta2*exp(-(tnow)*theta3)) 
      # no process noise for mu_t (which is i.i.d.)
      meer[p,i] = sqrt(R)*rnorm(1) 
      if (i > 1){
        e[p,i] = AR*e[p,i-1] + meer[p,i]
        # here is process noise (i.i.v - "independently and incrementally varied")
      }else{e[p,i] = meer[p,i]}
      y[p,i] = mu_t+ e[p,i]
    }
    y[p,(nT+1):(nT*2)] = time
    #  lines(time,y[p,1:nT],lwd=1,lty=1,col=1)
  }
  
  y = cbind(1:nP,y)
  #head(y)
  
  # reshape data from wide to long
  dat <- as.data.frame(y)
  colnames(dat) <- c("id", paste0("y",1:nT),paste0("time",1:nT))
  dat_long <- matrix(NA, nrow=nP*nT, ncol=(2+ydim))
  colnames(dat_long) <- c("id", "time", paste0("y",1:ydim))
  for(i in 1:nP){
    dat_long[((i-1)*nT+1):(i*nT), 1] = i
    for(t in 1:nT){
      dat_long[((i-1)*nT+t),2] = dat[i,1+nT+t]
      dat_long[((i-1)*nT+t),3] = dat[i,1+t]
    }
  }
  head(dat_long)
  
  #------------------- Save Data -------------------
  write.table(dat_long,paste0(dat_name, ".dat") ,row.names = F, col.names = F, quote = F)
  write.csv(dat_long,paste0(dat_name, ".csv"), row.names = F, col.names = T, quote = F)
  #write.table(dat,paste0(dat_name,"_wide", ".dat") ,row.names = F, col.names = F, quote = F)
  #write(t(y),file = paste0(dat_name,"_wide", ".dat"),ncol=dim(y)[2],append=FALSE)
  print(paste0(dat_name,".dat"))
  print(paste0(dat_name,".csv"))
  #print(paste0(dat_name,"_wide", ".dat"))
  
  #------------------- Check Data -------------------
  # # Comparison of trend-only data
  # y_trendonly = matrix(NA, nP, nT*2) #data matrix
  # e_trendonly = matrix(NA, nP, nT*2) #residual data matrix
  # par(mfrow=c(1,1))
  # plot(1:nT, 1:nT, type="n", xlim=c(0,11), 
  #      ylim=c(0,60), xlab="Time", ylab = "Observed trend-only data")
  # for (p in 1:nP){
  #   time = dat_long[((p-1)*nT+1):(p*nT),2]
  #   theta1 = theta_s[p,1]
  #   theta2 = theta_s[p,2]
  #   theta3 = theta_s[p,3]
  #   for (i in 1:nT){
  #     tnow = time[i]
  #     mu_t = theta1*exp(-theta2*exp(-(tnow)*theta3)) 
  #     #e_trendonly[p,i] = meer[p,i]
  #     y_trendonly[p,i] = mu_t #+ e_trendonly[p,i]
  #   }
  #   y_trendonly[p,(nT+1):(nT*2)] = t(time)
  #   lines(time,y_trendonly[p,1:nT],lwd=1,lty=1,col=1)
  # }
  # 
  # y_trendonly = cbind(1:nP,y_trendonly)
  # #head(y)
  # 
  # # Comparison of residual-only data
  # y_resonly = matrix(NA, nP, nT*2) #data matrix
  # par(mfrow=c(1,1))
  # plot(1:nT, 1:nT, type="n", xlim=c(0,11), 
  #      ylim=c(0,60), xlab="Time", ylab = "Observed residual-only data")
  # for (p in 1:nP){
  #   time = dat_long[((p-1)*nT+1):(p*nT),2]
  #   mu_t = mean(y_trendonly[p,1:nT+1])
  #   for (i in 1:nT){
  #     y_resonly[p,i] = mu_t+ e[p,i]
  #   }
  #   y_resonly[p,(nT+1):(nT*2)] = t(time)
  #   lines(time,y_resonly[p,1:nT],lwd=1,lty=1,col=1)
  # }
  # 
  # y_resonly = cbind(1:nP,y_resonly)
  
  # Comparison of measurement error tracess
  # par(mfrow=c(1,1))
  # plot(1:nT, 1:nT, type="n", xlim=c(0,11), 
  #      ylim=c(0,60), xlab="Time", ylab = "Observed measurement error")
  # for (p in 1:nP){
  #   time = dat_long[((p-1)*nT+1):(p*nT),2]
  #   mu_t = mean(y[p,1:nT+1])
  #   lines(time,mu_t+meer[p,1:nT],lwd=1,lty=1,col=1)
  # }
}

########################################################################
######################## MC Experiments ends ########################
########################################################################
Total_t_end <- proc.time()
Total_t_elapsed <- Total_t_end - Total_t_begin
show(Total_t_elapsed)
