rm(list=ls())
library(ggplot2)
library(ggsci) #https://zhuanlan.zhihu.com/p/430107483
library(cowplot)
library(dplyr)
theme_set(theme_gray(base_family="Arial Unicode MS"))

######################## Global Control Parameters ########################
work_path <- "~/XXYDATAanalysis/GoHiARmodel/substantivedata"
setwd(work_path)

C="ECLSK"
dateMark="0709"

sigmoid <- function(time, theta1, theta2, theta3){
  theta1*exp(-theta2*exp(-(time)*theta3)) - 3 ### so that the function rage: (-3, +infinite)
}

Mean = function(x){
  mean(x,na.rm=TRUE)}

SD = function(x){
  sd(x,na.rm=TRUE)}

source("posteriorSummaryStats.R")

condis = c("DSEMonly","Linear","NLSLIST","NLME") #singlestage?

summ_colnames = c("Est","SE")
summ_rownames = c("IIV",
                  "Level2MeanAR","Level2Var_AR",
                  "Level2Mean_theta1","Level2Mean_theta2","Level2Mean_theta3",
                  "Level2Var_theta1","Level2Var_theta2","Level2Var_theta3")
nT=5
######################## Read in Result Files | Perturbation ########################
outputs = data.frame()

# DSEMonly # check nP already
Condi = "DSEMonly"
nP=2369 # this is 2369 because the file name is generate at the beginging of script. the actual nP=2143
r=2023
alloutput_DSEMonly = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"))
#alloutput_DSEMonly$parameters$unstandardized
outputs[1,"condi"] = Condi
outputs[1,"nT"] = nT
outputs[1,"nP"] = 2143
outputs[1,"Est_IIV"] = alloutput_DSEMonly$parameters$unstandardized[1,"est"]
outputs[1,"Est_Level2MeanAR"]=alloutput_DSEMonly$parameters$unstandardized[3,"est"]
outputs[1,"Est_Level2Var_AR"]=alloutput_DSEMonly$parameters$unstandardized[5,"est"]
outputs[1,"se_IIV"] = alloutput_DSEMonly$parameters$unstandardized[1,"posterior_sd"]
outputs[1,"se_Level2MeanAR"]=alloutput_DSEMonly$parameters$unstandardized[3,"posterior_sd"]
outputs[1,"se_Level2Var_AR"]=alloutput_DSEMonly$parameters$unstandardized[5,"posterior_sd"]
outputs[1,"LL_IIV"] = alloutput_DSEMonly$parameters$unstandardized[1,"lower_2.5ci"]
outputs[1,"LL_Level2MeanAR"]=alloutput_DSEMonly$parameters$unstandardized[3,"lower_2.5ci"]
outputs[1,"LL_Level2Var_AR"]=alloutput_DSEMonly$parameters$unstandardized[5,"lower_2.5ci"]
outputs[1,"UL_IIV"] = alloutput_DSEMonly$parameters$unstandardized[1,"upper_2.5ci"]
outputs[1,"UL_Level2MeanAR"]=alloutput_DSEMonly$parameters$unstandardized[3,"upper_2.5ci"]
outputs[1,"UL_Level2Var_AR"]=alloutput_DSEMonly$parameters$unstandardized[5,"upper_2.5ci"]


# NLSLIST
Condi = "NLSLIST"
nP=2143
r=2023

alloutput_NLSLIST = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"))
outputs[2,"condi"] = Condi
outputs[2,"nT"] = nT
outputs[2,"nP"] = 2143
outputs[2,"Est_IIV"] = alloutput_NLSLIST$parameters$unstandardized[1,"est"]
outputs[2,"Est_Level2MeanAR"]=alloutput_NLSLIST$parameters$unstandardized[3,"est"]
outputs[2,"Est_Level2Var_AR"]=alloutput_NLSLIST$parameters$unstandardized[5,"est"]
outputs[2,"se_IIV"] = alloutput_NLSLIST$parameters$unstandardized[1,"posterior_sd"]
outputs[2,"se_Level2MeanAR"]=alloutput_NLSLIST$parameters$unstandardized[3,"posterior_sd"]
outputs[2,"se_Level2Var_AR"]=alloutput_NLSLIST$parameters$unstandardized[5,"posterior_sd"]
outputs[2,"LL_IIV"] = alloutput_NLSLIST$parameters$unstandardized[1,"lower_2.5ci"]
outputs[2,"LL_Level2MeanAR"]=alloutput_NLSLIST$parameters$unstandardized[3,"lower_2.5ci"]
outputs[2,"LL_Level2Var_AR"]=alloutput_NLSLIST$parameters$unstandardized[5,"lower_2.5ci"]
outputs[2,"UL_IIV"] = alloutput_NLSLIST$parameters$unstandardized[1,"upper_2.5ci"]
outputs[2,"UL_Level2MeanAR"]=alloutput_NLSLIST$parameters$unstandardized[3,"upper_2.5ci"]
outputs[2,"UL_Level2Var_AR"]=alloutput_NLSLIST$parameters$unstandardized[5,"upper_2.5ci"]

# NLME # check nP already
Condi = "NLME" # this is 2369 because the file name is generate at the beginning of script. the actual nP=2143
nP=2369
r=2023

alloutput_NLME = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"))
outputs[3,"condi"] = Condi
outputs[3,"nT"] = nT
outputs[3,"nP"] = 2369
outputs[3,"Est_IIV"] = alloutput_NLME$parameters$unstandardized[1,"est"]
outputs[3,"Est_Level2MeanAR"]=alloutput_NLME$parameters$unstandardized[3,"est"]
outputs[3,"Est_Level2Var_AR"]=alloutput_NLME$parameters$unstandardized[5,"est"]
outputs[3,"se_IIV"] = alloutput_NLME$parameters$unstandardized[1,"posterior_sd"]
outputs[3,"se_Level2MeanAR"]=alloutput_NLME$parameters$unstandardized[3,"posterior_sd"]
outputs[3,"se_Level2Var_AR"]=alloutput_NLME$parameters$unstandardized[5,"posterior_sd"]
outputs[3,"LL_IIV"] = alloutput_NLME$parameters$unstandardized[1,"lower_2.5ci"]
outputs[3,"LL_Level2MeanAR"]=alloutput_NLME$parameters$unstandardized[3,"lower_2.5ci"]
outputs[3,"LL_Level2Var_AR"]=alloutput_NLME$parameters$unstandardized[5,"lower_2.5ci"]
outputs[3,"UL_IIV"] = alloutput_NLME$parameters$unstandardized[1,"upper_2.5ci"]
outputs[3,"UL_Level2MeanAR"]=alloutput_NLME$parameters$unstandardized[3,"upper_2.5ci"]
outputs[3,"UL_Level2Var_AR"]=alloutput_NLME$parameters$unstandardized[5,"upper_2.5ci"]

# Linear # check nP already
Condi = "Linear" # this is 2369 because the file name is generate at the beginning of script. the actual nP=2143
nP=2369
r=2023

alloutput_Linear = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"))
outputs[4,"condi"] = Condi
outputs[4,"nT"] = nT
outputs[4,"nP"] = 2143
outputs[4,"Est_IIV"] = alloutput_Linear$parameters$unstandardized[1,"est"]
outputs[4,"Est_Level2MeanAR"]=alloutput_Linear$parameters$unstandardized[3,"est"]
outputs[4,"Est_Level2Var_AR"]=alloutput_Linear$parameters$unstandardized[5,"est"]
outputs[4,"se_IIV"] = alloutput_Linear$parameters$unstandardized[1,"posterior_sd"]
outputs[4,"se_Level2MeanAR"]=alloutput_Linear$parameters$unstandardized[3,"posterior_sd"]
outputs[4,"se_Level2Var_AR"]=alloutput_Linear$parameters$unstandardized[5,"posterior_sd"]
outputs[4,"LL_IIV"] = alloutput_Linear$parameters$unstandardized[1,"lower_2.5ci"]
outputs[4,"LL_Level2MeanAR"]=alloutput_Linear$parameters$unstandardized[3,"lower_2.5ci"]
outputs[4,"LL_Level2Var_AR"]=alloutput_Linear$parameters$unstandardized[5,"lower_2.5ci"]
outputs[4,"UL_IIV"] = alloutput_Linear$parameters$unstandardized[1,"upper_2.5ci"]
outputs[4,"UL_Level2MeanAR"]=alloutput_Linear$parameters$unstandardized[3,"upper_2.5ci"]
outputs[4,"UL_Level2Var_AR"]=alloutput_Linear$parameters$unstandardized[5,"upper_2.5ci"]

# Single-Stage
Condi = "Single_Stage_2_large"
nP=2369 # this is 2369 because the file name is generate at the beginning of script. the actual nP=2143
r=2023
file_level2pars = paste0(Condi,"_level2pars","_r",r,"_",dateMark,".csv")
jagsoutput = read_csv(file_level2pars)
EstParNames = paste0("Est_", summ_rownames)
seParNames = paste0("se_", summ_rownames)
LLParNames = paste0("LL_", summ_rownames)
ULParNames = paste0("UL_", summ_rownames)

outputs[5,"condi"] = Condi
outputs[5,EstParNames]=jagsoutput$mean
outputs[5,seParNames]=jagsoutput$PSD
outputs[5,LLParNames]=jagsoutput$`PCI 2.50%`
outputs[5,ULParNames]=jagsoutput$`PCI 97.50%`

######################## Read in Result Files | Trend ########################
# NLME
Condi = "NLME" # this is 2369 because the file name is generate at the beginning of script. the actual nP=2143
nP=2369
r=2023
nlmeoutput = readRDS(paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds"))
outputs[3,"Est_Level2Mean_theta1"] = nlmeoutput$coefficients$fixed[1]
outputs[3,"Est_Level2Mean_theta2"] = nlmeoutput$coefficients$fixed[2]  
outputs[3,"Est_Level2Mean_theta3"] = nlmeoutput$coefficients$fixed[3]

outputs[3,"Est_Level2Var_theta1"] = as.numeric(VarCorr(nlmeoutput)[1,1])  
outputs[3,"Est_Level2Var_theta2"] = as.numeric(VarCorr(nlmeoutput)[2,1])  
outputs[3,"Est_Level2Var_theta3"] = as.numeric(VarCorr(nlmeoutput)[3,1])

# NLSLIST
Condi = "NLSLIST"
nP=2143
r=2023
nlsListoutput = readRDS(paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds"))
outputs[2,"Est_Level2Mean_theta1"]=Mean(nlsListoutput$coefficients[,1,1])
outputs[2,"Est_Level2Mean_theta2"]=Mean(nlsListoutput$coefficients[,1,2])
outputs[2,"Est_Level2Mean_theta3"]=Mean(nlsListoutput$coefficients[,1,3])

outputs[2,c("Est_Level2Var_theta1")]=SD(nlsListoutput$coefficients[,1,1])^2
outputs[2,c("Est_Level2Var_theta2")]=SD(nlsListoutput$coefficients[,1,2])^2
outputs[2,c("Est_Level2Var_theta3")]=SD(nlsListoutput$coefficients[,1,3])^2


write.csv(outputs,"ECLS-K_summary.csv",col.names = T)

