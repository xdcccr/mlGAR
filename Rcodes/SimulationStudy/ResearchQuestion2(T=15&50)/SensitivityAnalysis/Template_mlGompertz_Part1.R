rm(list=ls())

######################## Set Environment ########################
# Define a list of required packages
required_packages <- c("mvtnorm", "Matrix", "dplyr", "rjags", 
                       "coda", "tidyr", "lubridate", "tidyverse",
                       "data.table", "nlme")

# Load all required packages
print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

work_path <- "D:/XXYDATAanalysis/IP_1/resultsRevision1/results_1"
print(paste0("work_path: ", work_path))
setwd(work_path)
source("posteriorSummaryStats.R")

######################## Control Parameters ########################
toPlot = FALSE 
dateMark = "0614"
C = "2_nonlinear"
Condi = "TwoStageNLME"

# Add sigmoid function definition
sigmoid <- function(time, theta1, theta2, theta3){
  theta1*exp(-theta2*exp(-(time)*theta3))
}

ydim = 1 # number of dimensions of within-level variables
nT = nT_REPLACE # number of observations in time series 
npad = 0  # initial number of time points to discard before retaining
nP = nP_REPLACE # number of people 

N_para = 9 # number of parameters
N_repl= 5 # number of replications | replication k 

Round = Round_REPLACE

rs = 1:N_repl+Round*N_repl

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

######################## Model Parameters ########################
#person-specific auto-regression coefficients
AR_mean <-  .3 # [mean of AR] 
AR_sigma <- .1 # [variances of AR] 
AR <- AR_mean + rnorm(nP,0,AR_sigma)

theta10 = 35 # asymptote
theta20 = 4  # displacement
theta30 = .8 # growth rate
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
for(k in 1:N_repl){ # k=1
  
  r=rs[k]
  set.seed(r)
  
  dat_name <- paste0("Data_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark)
  
  dat_filename = paste0(dat_name,".csv")
  file_detrdat_long = paste0(Condi,"_detrdat_long_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".txt")
  file_ROPEsamples = paste0(Condi,"_ROPEsamples_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  file_level2pars = paste0(Condi,"_level2pars_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  file_jagsoutput = paste0(Condi,"_jagsoutput_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  
  print(paste0("dat_filename: ", dat_filename))
  
  ######################## First Stage: Detrend Data ########################
  # Load data - long
  dat_fit <- read.csv(dat_filename)
  
  # Get dimensions
  ids = unique(dat_fit$id)
  
  # Detrend Data using NLME
  y1_sigmoid.nlme = nlme(y1 ~ sigmoid(time, theta1, theta2, theta3),
                         data = dat_fit,
                         fixed = theta1 + theta2 + theta3 ~ 1,
                         random = (theta1 + theta2 + theta3 ~ 1 | id),
                         start = c(theta1=theta10, 
                                   theta2=theta20, 
                                   theta3=theta30))
  
  res_detr <-summary(y1_sigmoid.nlme)
  fitted_nlme <- fitted(res_detr)
  resid_nlme  <- dat_fit$y1-fitted_nlme
  
  dat_detr <- matrix(NA, nrow=nrow(dat_fit), ncol=ncol(dat_fit))
  colnames(dat_detr) <- colnames(dat_fit)
  dat_detr[,c("id","time")]=as.matrix(dat_fit[,c("id","time")])
  dat_detr[,"y1"]=resid_nlme 
  dat_detr <- as.data.frame(dat_detr)
  
  # Save nlme results
  saveRDS(res_detr, file = paste0(Condi,"_nlme_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds"))
  write.table(dat_detr, file_detrdat_long, row.names = FALSE, col.names = TRUE, quote = F)
  
  ######################## Second Stage: JAGS Model ########################
  # Prepare data for JAGS
  # Convert long format to array format for JAGS
  Obs_empty = rep(NA, nP*nT*(ydim+1))
  Obs = array(Obs_empty,dim=c(nT, (ydim+1), nP))
  for (i in ids){ #i=2
    subset_i = subset(dat_detr,id==i)
    Obs[,,i] = as.matrix(subset_i[,-1])
  }
  
  # Create JAGS data list
  jagsData <- list("ids" = ids, "nT" = nT, "Obs" = Obs)
  
  # Specify the model for detrended data
  cat("
  model {
    # PART 1. Specifying Level-1 model
    for (id in ids) { # opening loop for person (id)
      # model for the first point
      Obs[1,2,id] ~ dnorm(0, 0.1)
      
      # model for the rest
      for (t in 2:nT) { 
        Obs[t,2,id] ~ dnorm(mu[id,t], 1/IIV)
        mu[id,t] <- AR[id]*(Obs[t-1,2,id])  
      }
      
      # Level-2 model
      AR[id] ~ dnorm(Level2MeanAR, 1/pow(Level2Sd_AR,2))T(-1,1)
    }
      
    # PART 2. Specifying prior distributions 
    # AR - using informative priors matching mlGAR
    Level2MeanAR ~ dnorm(0,1)T(-1,1)
    Level2Sd_AR ~ dunif(0,1)
    Level2Var_AR = pow(Level2Sd_AR,2)
    
    # IIV
    IIV <- exp(logIIV)   
    logIIV ~ dunif(-4.605170, 1.83) #log(c(0.1^2,2.5^2))
  }
  ", file = file_model)
  
  # Parameters to track
  parameters <- c("AR", "Level2MeanAR", "Level2Var_AR", "IIV")
  
  # MCMC settings
  adaptation <- 50000  
  burnin <- 10000
  chains <- 2
  thin <- 5 
  niter <- 50000
  
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
  write.csv(resulttable, file=file_ROPEsamples)
  
  # Save level 2 parameters specifically
  level2_params <- resulttable[c("IIV","Level2MeanAR","Level2Var_AR"),]
  write.csv(level2_params, file=file_level2pars)
  
  # Save full JAGS output
  sumjags <- summary(codaSamples)
  write.csv(sumjags$statistics, file=file_jagsoutput)
  
}

########################################################################
######################## MC Experiments ends ########################
########################################################################
Total_t_end <- proc.time()
Total_t_elapsed <- Total_t_end - Total_t_begin
show(Total_t_elapsed)