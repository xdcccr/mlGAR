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
N_repl= 100 # number of replications | replication k 

# Round = Round_REPLACE
# rs = 1:N_repl+Round*N_repl
rs = 1:N_repl

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
  file_model = paste0(Condi,"_nlme_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")
  
  # Read JAGS results
  jagsoutput = read_csv(file_level2pars)
  nlmeoutput = readRDS(file_model)
  
  MCfile[k,"nT"]=nT
  MCfile[k,"nP"]=nP
  MCfile[k,"repl" ]=k
  
  # Save true parameters
  MCfile[k,TruParNames]=c(R,
                          AR_mean,AR_sigma^2,
                          theta10,theta20,theta30,
                          sigma1^2,sigma2^2,sigma3^2)
  
  # Save AR parameters from JAGS
  MCfile[k,"Est_IIV"] = jagsoutput$mean[jagsoutput$`...1`=="IIV"]
  MCfile[k,"Est_Level2MeanAR"]=jagsoutput$mean[jagsoutput$`...1`=="Level2MeanAR"]
  MCfile[k,"Est_Level2Var_AR"]=jagsoutput$mean[jagsoutput$`...1`=="Level2Var_AR"]
  
  MCfile[k,"se_IIV"] = jagsoutput$PSD[jagsoutput$`...1`=="IIV"]
  MCfile[k,"se_Level2MeanAR"]=jagsoutput$PSD[jagsoutput$`...1`=="Level2MeanAR"]
  MCfile[k,"se_Level2Var_AR"]=jagsoutput$PSD[jagsoutput$`...1`=="Level2Var_AR"]
  
  MCfile[k,"LL_IIV"] = jagsoutput$`PCI 2.50%`[jagsoutput$`...1`=="IIV"]
  MCfile[k,"LL_Level2MeanAR"]=jagsoutput$`PCI 2.50%`[jagsoutput$`...1`=="Level2MeanAR"]
  MCfile[k,"LL_Level2Var_AR"]=jagsoutput$`PCI 2.50%`[jagsoutput$`...1`=="Level2Var_AR"]
  
  MCfile[k,"UL_IIV"] = jagsoutput$`PCI 97.50%`[jagsoutput$`...1`=="IIV"]
  MCfile[k,"UL_Level2MeanAR"]=jagsoutput$`PCI 97.50%`[jagsoutput$`...1`=="Level2MeanAR"]
  MCfile[k,"UL_Level2Var_AR"]=jagsoutput$`PCI 97.50%`[jagsoutput$`...1`=="Level2Var_AR"]
  # Save coverage flags
  MCfile[k,"cFlag_IIV"]=ifelse((MCfile[k,"Tru_IIV"]>=MCfile[k,"LL_IIV"]) && 
                                 (MCfile[k,"Tru_IIV"]<=MCfile[k,"UL_IIV"]),1,0)
  MCfile[k,"cFlag_Level2MeanAR"]=ifelse((MCfile[k,"Tru_Level2MeanAR"]>=MCfile[k,"LL_Level2MeanAR"]) && 
                                          (MCfile[k,"Tru_Level2MeanAR"]<=MCfile[k,"UL_Level2MeanAR"]),1,0)
  MCfile[k,"cFlag_Level2Var_AR"]=ifelse((MCfile[k,"Tru_Level2Var_AR"]>=MCfile[k,"LL_Level2Var_AR"]) && 
                                          (MCfile[k,"Tru_Level2Var_AR"]<=MCfile[k,"UL_Level2Var_AR"]),1,0)
  
  # Save NLME parameters
  MCfile[k,c("Est_Level2Mean_theta1",
             "Est_Level2Mean_theta2",
             "Est_Level2Mean_theta3")]=nlmeoutput$coefficients$fixed
  MCfile[k,c("se_Level2Mean_theta1",
             "se_Level2Mean_theta2",
             "se_Level2Mean_theta3")]=sqrt(diag(nlmeoutput$varFix))
  MCfile[k,c("LL_Level2Mean_theta1",
             "LL_Level2Mean_theta2",
             "LL_Level2Mean_theta3")]=nlmeoutput$coefficients$fixed-1.96*sqrt(diag(nlmeoutput$varFix))
  MCfile[k,c("UL_Level2Mean_theta1",
             "UL_Level2Mean_theta2",
             "UL_Level2Mean_theta3")]=nlmeoutput$coefficients$fixed+1.96*sqrt(diag(nlmeoutput$varFix))
  MCfile[k,"cFlag_Level2Mean_theta1"]=ifelse((MCfile[k,"Tru_Level2Mean_theta1"]>=MCfile[k,"LL_Level2Mean_theta1"]) && (MCfile[k,"Tru_Level2Mean_theta1"]<=MCfile[k,"UL_Level2Mean_theta1"]),1,0)
  MCfile[k,"cFlag_Level2Mean_theta2"]=ifelse((MCfile[k,"Tru_Level2Mean_theta2"]>=MCfile[k,"LL_Level2Mean_theta2"]) && (MCfile[k,"Tru_Level2Mean_theta2"]<=MCfile[k,"UL_Level2Mean_theta2"]),1,0)
  MCfile[k,"cFlag_Level2Mean_theta3"]=ifelse((MCfile[k,"Tru_Level2Mean_theta3"]>=MCfile[k,"LL_Level2Mean_theta3"]) && (MCfile[k,"Tru_Level2Mean_theta3"]<=MCfile[k,"UL_Level2Mean_theta3"]),1,0)
  
  MCfile[k,"Est_Level2Var_theta1"] = as.numeric(VarCorr(nlmeoutput)[1,1])  
  MCfile[k,"Est_Level2Var_theta2"] = as.numeric(VarCorr(nlmeoutput)[2,1])  
  MCfile[k,"Est_Level2Var_theta3"] = as.numeric(VarCorr(nlmeoutput)[3,1])
  
}

write.csv(MCfile, file=file_MCfile)

######################## Compute Summary Statistics ########################
MCfile = read.csv(file_MCfile)
remab=T

if(remab==T){
  cols_to_replace <- EstParNames
  
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
#------------------- Compute summary statistics -------------------  
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