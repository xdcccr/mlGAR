########################################################################
# The script is used to generate MCSumm files from MCfiles.
########################################################################
rm(list=ls())

required_packages <- c("mvtnorm", "Matrix", "dplyr", 
                       "tidyverse", "data.table")

work_path <- "D:/XXY/IP_1/result0614_PC/results"
print(paste0("work_path: ", work_path))
setwd(work_path)

#### Common Setting ####
sigmoid <- function(time, theta1, theta2, theta3){
  theta1*exp(-theta2*exp(-(time)*theta3))}
  
Mean = function(x){
  mean(x,na.rm=TRUE)}
  
SD = function(x){
  sd(x,na.rm=TRUE)}
  
relBias = function(x,truex){
  mean(((x-truex))/truex,na.rm=T)}
  
RMSE = function(x,truex){
  sqrt(Mean(x-truex)^2)}

R = 1 # measurement error variance
AR_mean <-  .3 # [mean of AR] 
AR_sigma <- .1 # [variances of AR] 

theta10 = 35 # asymptote, limitation when time gets close to +infinite
theta20 = 4 # sets the displacement along the x-axis (translates the graph to the left or right). 
theta30 = .8 #sets the growth rate (y scaling)

r12 = 0 #cor(theta_1,theta_2)
r13 = 0 #cor(theta_1,theta_3)
r23 = 0 #cor(theta_2,theta_3)
sigma1 = 9 #SD(theta_1)
sigma2 = .5  #SD(theta_2)
sigma3 = .10  #SD(theta_3)


for(McCondis in c("TwoStageNLSLIST","TwoStageNLME","TwoStageLinear")){
  for(MCnT in c(5,15,50)){
    for(MCnP in c(150,500)){
      ######################## Control Parameters ########################
      toPlot = FALSE #Added this flag to plot
      dateMark = "0614"
      C = "2_nonlinear"
      Condi = McCondis
      
      
      ydim = 1 # number of dimensions of within-level variables
      nT = MCnT # number of observations in time series 
      npad = 0  # initial number of time points to discard before retaining
      nP = MCnP # number of people 
      
      N_para = 9 # number of parameters
      N_repl= 100 # number of replications | replication k 
      #N_exp_condi = 6 # number of experiment conditions
      
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

      ######################## Compute Summary Statistics ########################
      MCfile = read.csv(file_MCfile)
      cols_to_replace <- EstParNames
      
      for(col in summ_rownames){
        Tcol=paste0("Tru_",col)
        Ecol=paste0("Est_",col)
        secol=paste0("se_",col)
        LLcol=paste0("LL_",col)
        ULcol=paste0("UL_",col)
        cFlagcol=paste0("cFlag_",col)
        for(row in 1:nrow(MCfile)){
          if(!is.na(MCfile[row, Ecol])){
            if(MCfile[row, Ecol] > max(3*MCfile[row, Tcol],10+MCfile[row, Tcol])){
              MCfile[row, Ecol] <- NA}
          if(is.na(MCfile[row, Ecol])){
            MCfile[row, secol] <- NA
            MCfile[row, LLcol] <- NA
            MCfile[row, ULcol] <- NA
            MCfile[row, cFlagcol] <- NA
          }
          }
        }
      }
      
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
    }
  }
}