rm(list=ls())

# List of packages
required_packages <- c("tidyr", "MplusAutomation", "texreg")

# Function to check and install missing packages
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE)) {
#     install.packages(pkg, repos='http://cran.rstudio.com/')
#   }
# }

# Check and install packages
print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

work_path <- "~/XXYDATAanalysis/GoHiARmodel/result0614_mac/results"
print(paste0("work_path: ", work_path))
setwd(work_path)
source("posteriorSummaryStats.R")

######################## Control Parameters ########################
C = "2_nonlinear"
Condi = "DSEMonly"

dateMark = "0614"
nT = 15 # number of observations in time series
nP = 150 # number of people
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
#deltaT <- 1 # [Time intervals] can be seconds/days/weeks/months/...
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
  print(paste0("dat_filename: ", dat_filename))
  dat_long=read.csv(dat_filename)
#------------------- Read in data -------------------  
  # # reshape data from wide to long
  # dat_wide <- matrix(NA,nrow=nP,ncol=1+(ydim+1)*nT)
  # colnames(dat_wide) <- c("id", paste0("y",1:nT),paste0("time",1:nT))
  # for(p in 1:nP){
  #   dat_wide[p, 1] = p
  #   dat_wide[p, 1+(1:nT)] = t(dat_long[((p-1)*nT+1):(p*nT),"y1"])
  #   dat_wide[p, 1+nT+(1:nT)] = t(dat_long[((p-1)*nT+1):(p*nT),"time"])
  # }
#------------------- Generate Mplus Script -------------------
  #Set Mplus input File
  input = mplusObject(TITLE = "A two-level DSEM model for one continuous dependent
                      variable with random intercepts and random slopes", 
                      VARIABLE = "MISSING ARE ALL (-99999);
                      CLUSTER = id;
                      LAGGED =  y1(1);",
                      ANALYSIS = "TYPE = TWOLEVEL RANDOM;
                      ESTIMATOR = BAYES; 
                      ALGORITHM = GIBBS;
                      PROCESSORS = 2;
                      BITERATIONS = 60000(10000);", 
                      MODEL = "%WITHIN%
                      s1 | y1 ON y1&1;

                      !%BETWEEN%
                      !y1 s1 WITH y1 s1;",
                      # MODELCONSTRAINT="NEW (s1_sd);
                      # s1_sd
                      # "
                      OUTPUT = "TECH1 TECH8;", #STANDARDIZED (CLUSTER)
                      SAVEDATA = paste0("SAVE=FSCORES(50 10);FILE IS ",Condi,"_FSCORES",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv;"),
                      PLOT = "TYPE = PLOT3;
                      FACTOR = ALL;",
                      usevariables = c("id","y1"),
                      rdata = dat_long,
                      autov = TRUE)
  res = mplusModeler(input, 
                     modelout = paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".inp"), 
                     hashfilename=FALSE,
                     run = 1L)
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
for(k in 1:N_repl){
  r=rs[k]
  allOutput = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"),recursive=TRUE)
  
  MCfile[k,"nT"]=nT
  MCfile[k,"nP"]=nP
  MCfile[k,"repl" ]=k
  
  MCfile[k,TruParNames]=c(R,
                          AR_mean,AR_sigma^2,
                          theta10,theta20,theta30,
                          sigma1^2,sigma2^2,sigma3^2)
  MCfile[k,"Est_IIV"] = allOutput$parameters$unstandardized[1,"est"]
  MCfile[k,"Est_Level2MeanAR"]=allOutput$parameters$unstandardized[3,"est"]
  MCfile[k,"Est_Level2Var_AR"]=allOutput$parameters$unstandardized[5,"est"]
  MCfile[k,"se_IIV"] = allOutput$parameters$unstandardized[1,"posterior_sd"]
  MCfile[k,"se_Level2MeanAR"]=allOutput$parameters$unstandardized[3,"posterior_sd"]
  MCfile[k,"se_Level2Var_AR"]=allOutput$parameters$unstandardized[5,"posterior_sd"]
  MCfile[k,"LL_IIV"] = allOutput$parameters$unstandardized[1,"lower_2.5ci"]
  MCfile[k,"LL_Level2MeanAR"]=allOutput$parameters$unstandardized[3,"lower_2.5ci"]
  MCfile[k,"LL_Level2Var_AR"]=allOutput$parameters$unstandardized[5,"lower_2.5ci"]
  MCfile[k,"UL_IIV"] = allOutput$parameters$unstandardized[1,"upper_2.5ci"]
  MCfile[k,"UL_Level2MeanAR"]=allOutput$parameters$unstandardized[3,"upper_2.5ci"]
  MCfile[k,"UL_Level2Var_AR"]=allOutput$parameters$unstandardized[5,"upper_2.5ci"]
  MCfile[k,"cFlag_IIV"]=ifelse((MCfile[k,"Tru_IIV"]>=MCfile[k,"LL_IIV"]) && (MCfile[k,"Tru_IIV"]<=MCfile[k,"UL_IIV"]),1,0)
  MCfile[k,"cFlag_Level2MeanAR"]=ifelse((MCfile[k,"Tru_Level2MeanAR"]>=MCfile[k,"LL_Level2MeanAR"]) && (MCfile[k,"Tru_Level2MeanAR"]<=MCfile[k,"UL_Level2MeanAR"]),1,0)
  MCfile[k,"cFlag_Level2Var_AR"]=ifelse((MCfile[k,"Tru_Level2Var_AR"]>=MCfile[k,"LL_Level2Var_AR"]) && (MCfile[k,"Tru_Level2Var_AR"]<=MCfile[k,"UL_Level2Var_AR"]),1,0)
}

write.csv(MCfile, file=file_MCfile)

#------------------- Compute summary statistics -------------------  
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
