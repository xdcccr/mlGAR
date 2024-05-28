########################################################################
# nlme: Nonlinear Mixed-Effects Models
  # This generic function fits a nonlinear mixed-effects model in 
    # the formulation described in Lindstrom and Bates (1990) 
    # but allowing for nested random effects. 
  # The within-group errors are allowed to be correlated and/or have unequal variances.
########################################################################

rm(list=ls())

# List of packages
required_packages <- c("mvtnorm", "Matrix", "dplyr", "tidyverse", 
                       "data.table", "forestmangr", "MplusAutomation", 
                       "texreg")

# Function to check and install missing packages
# for (pkg in required_packages) {
#   if (!require(pkg, character.only = TRUE)) {
#     install.packages(pkg, repos='http://cran.rstudio.com/')
#   }
# }

# Check and install packages
print("Load all required packages:")
lapply(required_packages, require, character.only = TRUE)

work_path <- "D:/XXYDATAanalysis/IP_1/result0614_PC/results"
print(paste0("work_path: ", work_path))
setwd(work_path)

######################## Control Parameters ########################
toPlot = FALSE #Added this flag to plot
dateMark = "0614"
C = "2_nonlinear"
Condi = "TwoStageLinear"

ydim = 1 # number of dimensions of within-level variables
nT = 5 # number of observations in time series 
npad = 0  # initial number of time points to discard before retaining
nP = 500 # number of people 

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

print(paste0("file_MCfile: ", file_MCfile))
print(paste0("file_MCfile_summ: ", file_MCfile_summ))

######################## Model Parameters ########################
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

theta10_start = 30 # asymptote, limitation when time gets close to +infinite
theta20_start = 2 # sets the displacement along the x-axis (translates the graph to the left or right).
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
  file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")
  # file_detrdat_wide = paste0(Condi,"_detrdat_wide_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv")
  

  print(paste0("dat_filename: ", dat_filename))
  print(paste0("file_detrdat_long: ", file_detrdat_long))
  print(paste0("file_model: ", file_model))
  # print(paste0("file_detrdat_wide: ", file_detrdat_long))

  ######################## Detrend Data ########################
  # Load data - long
  dat_fit <- read.csv(dat_filename)
  colnames(dat_fit)
  
  # Get dimensions
  # ids = unique(dat_fit$id)
  
  # Detrend Data
  y1_lm= lm_table(dat_fit, y1 ~ time, "id")
  
  dat_detr <- matrix(NA, nrow=nrow(dat_fit), ncol=ncol(dat_fit))
  colnames(dat_detr) <- colnames(dat_fit)
  dat_detr[,c("id","time")]=as.matrix(dat_fit[,c("id","time")])
  for(i in 1:nP){
    b0_i=y1_lm$b0[i]
    b1_i=y1_lm$b1[i]
    for(t in 1:nT){
      dat_detr[(i-1)*nT+t,"y1"] = dat_fit[(i-1)*nT+t,"y1"]-b0_i-b1_i*dat_detr[(i-1)*nT+t,"time"] 
    }
  }
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
  saveRDS(y1_lm, file = file_model)
  write.table(dat_detr_long, file_detrdat_long, row.names = FALSE, col.names = TRUE, quote = F)
  ######################## Fit DSEM in Mplus ########################
  # Set Mplus input File
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
                      %BETWEEN%
                      y1@0;
                      [y1@0];",
                      OUTPUT = "TECH1 TECH8;", #STANDARDIZED (CLUSTER)
                      SAVEDATA = paste0("SAVE=FSCORES(50 10);FILE IS ",Condi,"_FSCORES",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv;"),
                      PLOT = "TYPE = PLOT3;
                      FACTOR = ALL;",
                      usevariables = c("id","y1"),
                      rdata = dat_detr_long,
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

#------------------- Save necessary statistics -------------------  
print(dnames)
for(k in 1:N_repl){#k=1
  r=rs[k]
  file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")
  
  allOutput = readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"),recursive=TRUE)
  nlmeoutput = readRDS(file_model)
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

######################## Compute Summary Statistics ########################
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

#------------------- Compute summary statistics -------------------  
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