rm(list=ls())

######################## Set Environment ########################
# Define a list of required packages
required_packages <- c("mvtnorm", "Matrix", "dplyr", "rjags", 
                       "coda", "tidyr", "lubridate", "tidyverse",
                       "data.table", "nlme")

# Install any of the packages that are not yet installed
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE)) {
#     install.packages(pkg, repos='http://cran.rstudio.com/')
#   }
# }

# Load all required packages
print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

#work_path <- "/gpfs/group/quc16/default/xjx5093/IP_1/result0614"
work_path <- "~/XXYDATAanalysis/GoHiARmodel/result0614_mac/results"
print(paste0("work_path: ", work_path))
setwd(work_path)
source("posteriorSummaryStats.R")

######################## Control Parameters ########################
C = "2_nonlinear"
Condi = "SingleStagejags"

dateMark = "0614"
nT = 5 # number of observations in time series
nP = 500 # number of people
N_repl= 100 # number of replications | replication k
# N_exp_condi = 6 # number of experiment conditions

print(paste("dateMark:", dateMark, "nT:", nT, "nP:", nP, "N_repl:", N_repl))


ydim = 1 # number of dimensions of within-level variables

npad = 0  # initial number of time points to discard before retaining

N_para = 9 # number of parameters

rs = 1:N_repl

#------------------- MC files -------------------
TruParNames = paste0("TruPar", 1:N_para)
EstParNames = paste0("EstPar", 1:N_para)
seParNames = paste0("sePar", 1:N_para)

summ_colnames = c("TruePar","MeanPar_hat","RMSE","rBias",
                  "SD","Mean_SE_hat","RDSE",
                  "95%CI_LL","95%CI_UL","Coverage of 95% CIs")
summ_rownames = c("IIV",
                  "Level2MeanAR","Level2Var_AR",
                  "Level2Mean_theta1","Level2Mean_theta2","Level2Mean_theta3",
                  "Level2Var_theta1","Level2Var_theta2","Level2Var_theta3")

TruParNames = paste0("Tru_", summ_rownames)
EstParNames = paste0("Est_", summ_rownames)
seParNames = paste0("se_", summ_rownames)
LLParNames = paste0("LL_", summ_rownames)
ULParNames = paste0("UL_", summ_rownames)
CoverFlagNames = paste0("cFlag_", summ_rownames)

dnames = c("nT","nP","repl",
           TruParNames, EstParNames, seParNames,
           LLParNames,ULParNames,
           CoverFlagNames)

MCfile = matrix(NA,nrow = N_repl,ncol = length(dnames))
MCfile_summ = matrix(NA, nrow = N_para, ncol = length(summ_colnames)+1)

colnames(MCfile) = dnames
colnames(MCfile_summ) = c(summ_colnames,"MissingPer")
rownames(MCfile_summ) = summ_rownames

file_MCfile = paste0(Condi,"_MCfile_nT",nT,"_nP",nP,"_Nrepl",N_repl,"_",dateMark,".csv")
file_MCfile_summ = paste0(Condi,"_MCfileSumm_nT",nT,"_nP",nP,"_Nrepl",N_repl,"_",dateMark,".csv")

print(paste0("file_MCfile: ", file_MCfile))
print(paste0("file_MCfile_summ: ", file_MCfile_summ))

######################## Model Parameters ########################
# deltaT <- 1 # [Time intervals] can be seconds/days/weeks/months/...
sampleTime = seq(0.1,10,0.1)
#time <- matrix(seq(from = 0, to = (nT+npad-1)*deltaT, by = deltaT),
#               nrow = (nT + npad), ncol = 1) # regularly spaced

#person-specific auto-regression coefficients
AR_mean <-  .3 # [mean of AR] 
AR_sigma <- .1 # [variances of AR] 
AR <- AR_mean + rnorm(nP,0,AR_sigma)

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
for(k in 1:N_repl){

r=rs[k]
set.seed(r)
  
dat_name <- paste0("Data_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark)

dat_filename = paste0(dat_name,".csv")
file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".txt")
file_ROPEsamples = paste0(Condi,"_ROPEsamples_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
file_level2pars = paste0(Condi,"_level2pars_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
file_jagsoutput = paste0(Condi,"_jagsoutput_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")

print(paste0("dat_filename: ", dat_filename))
print(paste0("file_model: ", file_model))
print(paste0("file_ROPEsamples: ", file_ROPEsamples))
print(paste0("file_level2pars: ", file_level2pars))
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
parameters <- c("AR","MU","MU_thetas",
                summ_rownames)

# Specifying sampler settings
adaptation  <- 50000
burnin  <- 10000
chains  <- 2
thinning  <- 5
# Defining the number of posterior samples per chain for JAGS
nrOfIter <- 2000

# fixing the random seed for reproducibility
fixedinits<- list(list(.RNG.seed=6,.RNG.name="base::Mersenne-Twister"),
                  list(.RNG.seed=11,.RNG.name="base::Mersenne-Twister"))
main_t_start <- proc.time()

# Creating JAGS model object
jagsModel<-jags.model(file_model,data=jagsData,n.chains=chains,n.adapt=adaptation,inits=fixedinits)
# Running burn-in iterations
update(jagsModel,n.iter=burnin)
# Drawing posterior samples
codaSamples<-coda.samples(jagsModel,variable.names=parameters,n.iter=nrOfIter, thin = thinning,seed=5)

main_t_end <- proc.time()
main_t_elapsed <- main_t_end - main_t_start
show(main_t_elapsed)

# # Check specific sample - if needed
# abnorm_id = dat_fit[which(dat_fit$id==66),]
# y_limits <- range(dat_fit$y1)
# x_limits <- range(dat_fit$time)
# x_breaks <- seq(x_limits[1], x_limits[2], by = 1)
# ggplot(abnorm_id , aes(x = time, y = y1)) +
#   geom_point(size=0.5, color ="#00BFC4") + 
#   geom_line(color ="#00BFC4") +
#   scale_x_continuous(limits = x_limits, breaks = x_breaks) +
#   scale_y_continuous(limits = y_limits) +
#   ggtitle("abnormal individual trace") +
#   xlab("Time") +
#   ylab("y1") + theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     strip.background = element_blank(),
#     legend.position = "none"
#   )

######################## Checking Results ########################

parameters_2 <- summ_rownames

#------------------- Save statistics for selected parameters -------------------
resulttable <- summarizePost(codaSamples)
result_level_2 = resulttable[parameters_2,]
sumjags = summary(codaSamples)

write.csv(resulttable, file=file_ROPEsamples)
write.csv(result_level_2, file=file_level2pars)
write.csv(sumjags$statistics, file=file_jagsoutput)
}
########################################################################
######################## MC Experiments ends ########################
########################################################################
Total_t_end <- proc.time()
Total_t_elapsed <- Total_t_end - Total_t_begin
show(Total_t_elapsed)

######################## Save MC results ########################
#------------------- Build selection function -------------------  
Mean = function(x){
  mean(x,na.rm=TRUE)}

SD = function(x){
  sd(x,na.rm=TRUE)}

relBias = function(x,truex){
  mean(((x-truex))/truex,na.rm=T)
}

RMSE = function(x,truex){
  sqrt(Mean(x-truex)^2)
}

#------------------- Save necessary statistics -------------------  
print(dnames)
for(k in 1:N_repl){ #k=1
  r=rs[k]
  file_level2pars = paste0(Condi,"_level2pars_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  jagsoutput = read_csv(file_level2pars)
  
  
  MCfile[k,"nT"]=nT
  MCfile[k,"nP"]=nP
  MCfile[k,"repl" ]=k
  
  MCfile[k,TruParNames]=c(R,
                          AR_mean,AR_sigma^2,
                          theta10,theta20,theta30,
                          sigma1^2,sigma2^2,sigma3^2)
  MCfile[k,EstParNames]=jagsoutput$mean
  MCfile[k,seParNames]=jagsoutput$PSD
  MCfile[k,LLParNames]=jagsoutput$`PCI 2.50%`
  MCfile[k,ULParNames]=jagsoutput$`PCI 97.50%`
  for(par in summ_rownames){
    MCfile[k,paste0("cFlag_", par)]=
      ifelse((MCfile[k,paste0("Tru_", par)]>=MCfile[k,paste0("LL_", par)])
             && (MCfile[k,paste0("Tru_", par)]<=MCfile[k,paste0("UL_", par)]),1,0)
  }
}

write.csv(MCfile, file=file_MCfile)

#------------------- Compute Summary Statistics -------------------
MCfile = read.csv(file_MCfile)
remab=T

if(remab==T){
  for(col in summ_rownames){
    Tcol=paste0("Tru_",col)
    Ecol=paste0("Est_",col)
    for(row in 1:nrow(MCfile)){
      if(!is.na(MCfile[row, Ecol])){
        if(MCfile[row, Ecol] > max(3*MCfile[row, Tcol],10+MCfile[row, Tcol])){
          MCfile[row, Ecol] <- NA}
      }
    }
  }
}

print(summ_colnames)
print(summ_rownames)

MCfile_summ[,"TruePar"]=c(R,
                          AR_mean,AR_sigma^2,
                          theta10,theta20,theta30,
                          sigma1^2,sigma2^2,sigma3^2)
for(par in summ_rownames){#par="IIV"
  MCfile_summ[par,"MeanPar_hat"] = Mean(MCfile[,paste0("Est_",par)])
  MCfile_summ[par,"RMSE"] = RMSE(MCfile[,paste0("Est_",par)],
                                 MCfile[,paste0("Tru_",par)])
  MCfile_summ[par,"rBias"] = relBias(MCfile[,paste0("Est_",par)],
                                     MCfile[,paste0("Tru_",par)])
  MCfile_summ[par,"SD"] = SD(MCfile[,paste0("Est_",par)])
  MCfile_summ[par,"Mean_SE_hat"] = Mean(MCfile[,paste0("se_",par)])
  MCfile_summ[par,"RDSE"] = (MCfile_summ[par,"Mean_SE_hat"]-
                               MCfile_summ[par,"SD"])/MCfile_summ[par,"SD"]
  MCfile_summ[par,"RDSE"] = (MCfile_summ[par,"Mean_SE_hat"]-
                               MCfile_summ[par,"SD"])/MCfile_summ[par,"SD"]
  MCfile_summ[par,"95%CI_LL"]=quantile(MCfile[,paste0("Est_",par)], probs = 0.025, na.rm = T)
  MCfile_summ[par,"95%CI_UL"]=quantile(MCfile[,paste0("Est_",par)], probs = 0.975, na.rm = T)
  #MCfile_summ[par,"Coverage of 95% CIs"]= sum(MCfile[,paste0("cFlag_",par)],na.rm = T)/(N_repl-sum(is.na(MCfile[,paste0("cFlag_",par)])))
  MCfile_summ[par,"Coverage of 95% CIs"]= sum(MCfile[,paste0("cFlag_",par)],na.rm = T)/N_repl
  MCfile_summ[par,"MissingPer"]=sum(is.na(MCfile[,paste0("Est_",par)]))/(N_repl)
  }

write.csv(MCfile_summ, file=file_MCfile_summ)