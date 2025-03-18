rm(list=ls())
required_packages <- c("mvtnorm", "Matrix", "dplyr", 
                       "tidyverse", "data.table", "nlme", 
                       "rjags", "texreg")
# Check and install packages
print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

# work_path <- "/Users/xiaoyuexiong/XXYDATAanalysis/IP1-DetrendingPanelTimeSeries/CodesAndResultsForRevision/Everything2502"
work_path <- "D:/XXYDATAanalysis/IP_1/ClusterFile_PC/"
source("D:/XXYDATAanalysis/IP_1/ClusterFile_PC/posteriorSummaryStats.R")
print(paste0("work_path: ", work_path))
setwd(work_path)

######################## Control Parameters ########################
toPlot = FALSE #Added this flag to plot
dateMark = "2502"
C = "dim1_mlGAR_v2"
Condi = "mlGompertz"

ydim = 1 # number of dimensions of within-level variables
nT <- nT_REPLACE    # 第一个参数是nT
nP <- nP_REPLACE    # 第二个参数是nP
r <- r_REPLACE     # 第三个参数是replication number

# summ_colnames = c("TruePar","MeanPar_hat","RMSE","rBias",
#                   "SD","Mean_SE_hat","RDSE",
#                   "95%CI_LL","95%CI_UL","Coverage of 95% CIs")
summ_rownames = c("IIV",
                  "Level2MeanAR","Level2Var_AR",
                  "Level2Mean_theta1","Level2Mean_theta2","Level2Mean_theta3",
                  "Level2Var_theta1","Level2Var_theta2","Level2Var_theta3")

# TruParNames = paste0("Tru_", summ_rownames)
# EstParNames = paste0("Est_", summ_rownames)
# seParNames = paste0("se_", summ_rownames)
# LLParNames = paste0("LL_", summ_rownames)
# ULParNames = paste0("UL_", summ_rownames)
# CoverFlagNames = paste0("cFlag_", summ_rownames)

# dnames = c("nT","nP","repl",
#            TruParNames, EstParNames, seParNames,
#            LLParNames,ULParNames,
#            CoverFlagNames)
# 
# MCfile = matrix(NA,nrow = N_repl,ncol = length(dnames))
# MCfile_summ = matrix(NA, nrow = N_para, ncol = length(summ_colnames)+1)
# 
# colnames(MCfile) = dnames
# colnames(MCfile_summ) = c(summ_colnames,"MissingPer")
# rownames(MCfile_summ) = summ_rownames
# 
# file_MCfile = paste0(Condi,"_MCfile_nT",nT,"_nP",nP,"_Nrepl",N_repl,"_",dateMark,".csv")
# file_MCfile_summ = paste0(Condi,"_MCfileSumm_nT",nT,"_nP",nP,"_Nrepl",N_repl,"_",dateMark,".csv")
# 
# print(paste0("file_MCfile: ", file_MCfile))
# print(paste0("file_MCfile_summ: ", file_MCfile_summ))

######################## Model Parameters ########################
sampleTime = seq(0.1,10,0.1)

#person-specific auto-regression coefficients
AR_mean <-  .3 # [mean of AR] 
AR_sigma <- .1 # [variances of AR] 

y = matrix(NA, nP, nT*2) #data matrix
theta_s = matrix(NA, nP, 3) #data matrix
e = matrix(NA, nP, nT*2) #residual data matrix
meer = matrix(NA, nP, nT) #data matrix

theta10_start = 35 # asymptote, limitation when time gets close to +infinite
theta20_start = 4 # sets the displacement along the x-axis (translates the graph to the left or right).
theta30_start = 1 #sets the growth rate (y scaling)

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
sigmoid <- function(time, theta1, theta2, theta3){
  theta1*exp(-theta2*exp(-(time)*theta3))
}
######################## mlGompertz ########################
set.seed(r)

dat_name <- paste0("data/Data_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark)

dat_filename = paste0(dat_name,".csv")
file_detrdat_long = paste0(Condi,"_detrdat_long_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")

print(paste0("dat_filename: ", dat_filename))
print(paste0("file_detrdat_long: ", file_detrdat_long))
print(paste0("file_model: ", file_model))
# print(paste0("file_detrdat_wide: ", file_detrdat_long))

######################## Detrend Data ########################
# Load data - long
dat_fit <- read.csv(dat_filename)
colnames(dat_fit)

# Get dimensions
ids = unique(dat_fit$id)

# Detrend Data
y1_sigmoid.nlme = nlme(y1 ~ sigmoid(time, theta1, theta2, theta3),
                       data = dat_fit,
                       fixed = theta1 + theta2 + theta3 ~ 1,
                       random = (theta1 + theta2 + theta3 ~ 1 | id),
                       start = c(theta1=theta10_start, 
                                 theta2=theta20_start, 
                                 theta3=theta30_start))

# correlation = corAR1(AR_mean, form = ~1|id, fixed=FALSE)
res_detr <-summary(y1_sigmoid.nlme)
fitted_nlme <- fitted(res_detr)
resid_nlme  <- dat_fit$y1-fitted_nlme

dat_detr <- matrix(NA, nrow=nrow(dat_fit), ncol=ncol(dat_fit))
colnames(dat_detr) <- colnames(dat_fit)
dat_detr[,c("id","time")]=as.matrix(dat_fit[,c("id","time")])
dat_detr[,"y1"]=resid_nlme 
dat_detr <- as.data.frame(dat_detr)

if (toPlot==TRUE){
  # Check Data
  IDtoPlot = sort(sample(1:nP, 20, replace = FALSE))
  data_subset_1 = dat_fit[dat_fit$id %in% IDtoPlot,]
  data_subset_1$group <- "trended"
  data_subset_1_detrended =dat_detr[dat_detr$id %in% IDtoPlot,]
  data_subset_1_detrended$group="detrended"
  
  y_limits <- range(data_subset_1$y1)
  x_limits <- range(data_subset_1$time)
  x_breaks <- seq(x_limits[1], x_limits[2], by = 1)
  
  ggplot(data_subset_1, aes(x = time, y = y1, group = id)) +
    geom_point(size=0.5, color ="#00BFC4") +
    geom_line(color ="#00BFC4") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    scale_y_continuous(limits = y_limits) +
    ggtitle("individual traces") +
    xlab("Time") +
    ylab("y1") + theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      legend.position = "none"
    )
  
  ggplot(data_subset_1_detrended, aes(x = time, y = y1, group = id)) +
    geom_point(size=0.5, color="#F8766D") +
    geom_line(color="#F8766D") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    #scale_y_continuous(limits = y_limits) +
    ggtitle("individual traces") +
    xlab("Time") +
    ylab("y1") + theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      legend.position = "none"
    )
}

dat_detr_long <- as.data.frame(dat_detr)

# Save data 
saveRDS(res_detr, file = file_model)
write.table(dat_detr_long, file_detrdat_long, row.names = FALSE, col.names = TRUE, quote = F)


data_JAGS_equtime =  dat_detr_long
file_model = paste0(C,"_",Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".txt")
file_ROPEsamples = paste0(C,"_",Condi,"_ROPEsamples_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
# file_level2pars = paste0(C,"_",Condi,"_level2pars_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
file_jagsoutput = paste0(C,"_",Condi,"_jagsoutput_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")

nP=length(unique(data_JAGS_equtime$id))
nT=length(unique(data_JAGS_equtime$time))
ydim = 1
# Prepare data for JAGS
# Convert long format to array format for JAGS
Obs_empty = rep(NA, nP*nT*(ydim+1))
Obs = array(Obs_empty,dim=c(nT, (ydim+1), nP))
is_observed = array(Obs_empty,dim=c(nT, (ydim+1), nP))

for (p in 1:nP){ #i=2
  id_now = ids[p]
  subset_i = subset(data_JAGS_equtime,id==id_now)
  Obs[,,p] <- as.matrix(subset_i[,-1])
}

for(p in 1:nP){
  for(t in 1:nT){
    is_observed[t,2,p] = ifelse(!is.na(Obs[t,2,p]), 1, 0)
  }
}

# Create JAGS data list
jagsData <- list("ids" = ids, "nT" = nT, "nP"=nP,"Obs" = Obs,
                 "is_observed" = is_observed)
load.module("dic")
# Specify the model for detrended data
cat("
  model {
    # PART 1. Specifying Level-1 model
    for (p in 1:nP) { # opening loop for person (id)
      # model for the first point
      Obs[1,2,p] ~ dnorm(0, 0.1)
      log_lik[p,1] <- logdensity.norm(Obs[1,2,p], 0, 0.1)* is_observed[1,2,p]
      
      # model for the rest
      for (t in 2:nT) { 
        Obs[t,2,p] ~ dnorm(mu[p,t], 1/IIV)
        mu[p,t] <- AR[p]*(Obs[t-1,2,p])
        log_lik[p,t] <- logdensity.norm(Obs[t,2,p], mu[p,t], 1/IIV)* is_observed[t,2,p]
      }
      
      # Level-2 model
      AR[p] ~ dnorm(Level2MeanAR, 1/pow(Level2Sd_AR,2))T(-1,1)
    }
      
    # PART 2. Specifying prior distributions 
    # AR - using informative priors matching mlGAR
    Level2MeanAR ~ dnorm(0,1)T(-1,1)
    Level2Sd_AR ~ dunif(0,1)
    Level2Var_AR = pow(Level2Sd_AR,2)
    
    # IIV
    IIV <- exp(logIIV)   
    logIIV ~ dunif(-4.605170, 1.83) #log(c(0.1^2,2.5^2))
    
    # compute deviance
    total_log_lik <- sum(log_lik[1:nP, 1:nT])
    Dev <- -2 * total_log_lik
  }
  ", file = file_model)

# Parameters to track
parameters <- c("Level2MeanAR", "Level2Var_AR", "IIV", "Dev","deviance")

# MCMC settings
adaptation <- 5000 
burnin <- 5000
chains <- 2
thin <- 1 
niter <- 20000

# Set initial values
fixedinits <- list(
  list(.RNG.seed=6,.RNG.name="base::Mersenne-Twister"),
  list(.RNG.seed=11,.RNG.name="base::Mersenne-Twister")
)

# Fit JAGS model
jagsModel <- jags.model(file_model, 
                        data=jagsData, 
                        n.chains=chains,
                        n.adapt=adaptation,
                        inits=fixedinits)

update(jagsModel, n.iter=burnin)

codaSamples <- coda.samples(jagsModel,
                            variable.names=parameters,
                            n.iter=niter, 
                            thin=thin)

# Save results

resulttable <- summarizePost(codaSamples)
N_para=3
Estdeviance <- resulttable$mean[5]
Estbic<-Estdeviance+N_para*log(nT*nP)
resulttable <- rbind(resulttable,rep(NA,ncol(resulttable)))
rownames(resulttable)[6] <-"bic"
resulttable[6,1]=Estbic

write.csv(resulttable, file=file_ROPEsamples)

# # Save level 2 parameters specifically
# level2_params <- resulttable[c("IIV","Level2MeanAR","Level2Var_AR","Dev","deviance"),]
# write.csv(level2_params, file=file_level2pars)

# Save full JAGS output
sumjags <- summary(codaSamples)
write.csv(sumjags$statistics, file=file_jagsoutput)
