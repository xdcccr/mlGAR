rm(list=ls())

######################## Set Environment ########################
required_packages <- c("mvtnorm", "Matrix", "dplyr", "coda", "rjags", 
                       "tidyr", "lubridate", "tidyverse",
                       "data.table", "nlme")

print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

# 设置工作目录为相对路径
work_path <- "D:/XXYDATAanalysis/IP_1/ClusterFile_PC/"
print(paste0("work_path: ", work_path))
setwd(work_path)
source("D:/XXYDATAanalysis/IP_1/ClusterFile_PC/posteriorSummaryStats.R")
######################## 从命令行获取参数 ########################
nT <- nT_REPLACE    # 第一个参数是nT
nP <- nP_REPLACE    # 第二个参数是nP
r <- r_REPLACE     # 第三个参数是replication number

######################## Control Parameters ########################
C = "dim1_mlGAR_v2"
Condi = "SingleStageNoCor"
dateMark = "2502"
ydim = 1

print(paste("Processing: nT =", nT, "nP =", nP, "r =", r))

# #------------------- MC files -------------------
summ_rownames = c("IIV",
                  "Level2MeanAR","Level2Var_AR",
                  "Level2Mean_theta1","Level2Mean_theta2","Level2Mean_theta3",
                  "Level2Var_theta1","Level2Var_theta2","Level2Var_theta3")
N_para = length(summ_rownames) # number of parameters

######################## Model Parameters ########################
# deltaT <- 1 # [Time intervals] can be seconds/days/weeks/months/...
sampleTime = seq(0.1,10,0.1)
#time <- matrix(seq(from = 0, to = (nT+npad-1)*deltaT, by = deltaT),
#               nrow = (nT + npad), ncol = 1) # regularly spaced

#person-specific auto-regression coefficients
AR_mean <-  .3 # [mean of AR] 
AR_sigma <- .1 # [variances of AR] 
#AR <- AR_mean + rnorm(nP,0,AR_sigma)

y = matrix(NA, nP, nT*2) #data matrix
theta_s = matrix(NA, nP, 3) #data matrix
e = matrix(NA, nP, nT*2) #residual data matrix
meer = matrix(NA, nP, nT) #data matrix

theta10 = 35 # asymptote, limitation when time gets close to +infinite
theta20 = 4 # sets the displacement along the x-axis (translates the graph to the left or right). 
theta30 = .8 #sets the growth rate (y scaling)
r12 = 0 #cor(theta_1,theta_2)
r13 = 0 #cor(theta_1,theta_3)
r23 = 0 #cor(theta_2,theta_3)
sigma1 = 9 #SD(theta_1)
sigma2 = .5  #SD(theta_2)
sigma3 = .10  #SD(theta_3)
Q = matrix(c(sigma1^2, r12*sigma1*sigma2, r13*sigma1*sigma3,
             r12*sigma1*sigma2, sigma2^2, r23*sigma2*sigma3,
             r13*sigma1*sigma3, r23*sigma2*sigma3, sigma3^2), byrow=T,ncol=3)
Qr = chol(Q)

R = 1 # measurement error variance

########################################################################
######################## MC Experiments begins ########################
########################################################################
Total_t_begin <- proc.time()
  set.seed(r)
  
  dat_name <- paste0("data/Data_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark)
  
  dat_filename = paste0(dat_name,".csv")
  file_model = paste0("results/",Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".txt")
  file_ROPEsamples = paste0("results/",Condi,"_ROPEsamples_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  #file_level2pars = paste0("results/",Condi,"_level2pars_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  file_jagsoutput = paste0("results/",Condi,"_jagsoutput_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  
  print(paste0("dat_filename: ", dat_filename))
  print(paste0("file_model: ", file_model))
  print(paste0("file_ROPEsamples: ", file_ROPEsamples))
  #print(paste0("file_level2pars: ", file_level2pars))
  print(paste0("file_jagsoutput: ", file_jagsoutput))
  
  ######################## Model Fitting ########################
  # setwd(work_path)
  #------------------- Read and pre-process data -------------------
  # Load data - long
  dat_fit <- read.csv(dat_filename)
  colnames(dat_fit)
  
  # Get dimensions
  ids = unique(dat_fit$id)
  
  # Mutate the structure of dat to Y[nT,dim(y+time),id]
  Obs_empty = rep(NA, nP*nT*(ydim+1))
  Obs = array(Obs_empty,dim=c(nT, (ydim+1), nP))
  for (i in ids){ #i=2
    subset_i = subset(dat_fit,id==i)
    Obs[,,i] = as.matrix(subset_i[,-1])
  }
  Obs[,,1]  # data for id=1; columns = (time, y1)
  
  #------------------- MCMC -------------------
  # Create a list of all the variables 
  jagsData <- list("ids" = ids, "nT" = nT, "Obs" = Obs)
  load.module("dic")
  # Specify the Structured Latent Curve Model
  SLCmodel_0 = cat("
model {
# PART 1. Specifying Level-1 model
  for (id in ids) { # opening loop for person (id)
    # model for the first point (t=1)
    Obs[1,2,id] ~ dnorm(MU[id,1], 1/IIV)
    MU[id,1] <- MU_thetas[id,1]*exp(-MU_thetas[id,2]*exp(-Obs[1,1,id]*MU_thetas[id,3]))
    # model for the rest
    for (t in 2:nT) { # opening loop for observation (t)
      Obs[t,2,id] ~ dnorm(MU[id,t] + AR[id]*(Obs[t-1,2,id]-MU[id,t-1]), 1/IIV)
      MU[id,t] <- MU_thetas[id,1]*exp(-MU_thetas[id,2]*exp(-Obs[t,1,id]*MU_thetas[id,3]))
    } # closing loop for observation (t)
    
# PART 2. Specifying Level-2 model # should match data generation
  AR[id] ~ dnorm(Level2MeanAR, 1/pow(Level2Sd_AR,2))
  MU_thetas[id,1] ~ dnorm(Level2Mean_theta1, 1/pow(Level2Sd_theta1,2))
  MU_thetas[id,2] ~ dnorm(Level2Mean_theta2, 1/pow(Level2Sd_theta2,2)) 
  MU_thetas[id,3] ~ dnorm(Level2Mean_theta3, 1/pow(Level2Sd_theta3,2)) 
    
 } # closing loop for  person (id)       
      
# PART 3. Specifying prior distributions 
  # AR
    Level2MeanAR ~ dnorm(0,1)T(-1,1) #AR_mean = .2 # [Mean of AR] 
    Level2Sd_AR ~ dunif(0,1)  #AR_sigma = .1 # [variances of AR] 
    Level2Var_AR = pow(Level2Sd_AR,2)
  # MU
    Level2Mean_theta1 ~ dunif(20, 50) #theta10 = 35
    Level2Mean_theta2 ~ dunif(0, 15) #theta20 = 4
    Level2Mean_theta3 ~ dunif(0.5, 1.5) #theta30 = .8
    Level2Sd_theta1 ~ dunif(0.01,15) #sigma1 =  
    Level2Sd_theta2 ~ dunif(0.01,2) #sigma2 =   
    Level2Sd_theta3 ~ dunif(0.01,1) #sigma3 = 
    Level2Var_theta1 = pow(Level2Sd_theta1,2)
    Level2Var_theta2 = pow(Level2Sd_theta2,2)
    Level2Var_theta3 = pow(Level2Sd_theta3,2)
    
  # IIV # IIV = 1
    IIV <- exp(logIIV)   
    logIIV ~ dunif(-4.605170, 1.83) #log(c(0.1^2,2.5^2))
}
",file = file_model)
  
  # parameters to track
  parameters <- c(summ_rownames,"deviance")
  
  # Specifying sampler settings
  adaptation  <- 5000
  burnin  <- 5000
  chains  <- 2
  thinning  <- 1
  nrOfIter <- 20000
  
  # fixing the random seed for reproducibility
  fixedinits<- list(list(.RNG.seed=6,.RNG.name="base::Mersenne-Twister"),
                    list(.RNG.seed=11,.RNG.name="base::Mersenne-Twister"))
  main_t_start <- proc.time()
  
  # Creating JAGS model object
  jagsModel<-jags.model(file_model,
                        data=jagsData,
                        n.chains=chains,
                        n.adapt=adaptation,
                        inits=fixedinits)
  # Running burn-in iterations
  update(jagsModel,n.iter=burnin)
  # Drawing posterior samples
  codaSamples<-coda.samples(jagsModel,
                            variable.names=parameters,
                            n.iter=nrOfIter, 
                            thin = thinning,seed=5)
  
  main_t_end <- proc.time()
  main_t_elapsed <- main_t_end - main_t_start
  show(main_t_elapsed)
  ######################## Checking Results ########################
  
  #parameters_2 <- summ_rownames
  
  #------------------- Save statistics for selected parameters -------------------
  resulttable <- summarizePost(codaSamples)
  Estdeviance <- resulttable$mean[N_para+1]
  Estbic<-Estdeviance+N_para*log(nT*nP)
  resulttable <- rbind(resulttable,rep(NA,ncol(resulttable)))
  rownames(resulttable)[N_para+2] <-"bic"
  resulttable[N_para+2,1]=Estbic
  #result_level_2 = resulttable[parameters_2,]
  sumjags = summary(codaSamples)
  
  write.csv(resulttable, file=file_ROPEsamples)
  #write.csv(result_level_2, file=file_level2pars)
  write.csv(sumjags$statistics, file=file_jagsoutput)
########################################################################
######################## MC Experiments ends ########################
########################################################################
Total_t_end <- proc.time()
Total_t_elapsed <- Total_t_end - Total_t_begin
show(Total_t_elapsed)

