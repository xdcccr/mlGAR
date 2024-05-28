rm(list=ls())

# Data Cleaning
library(tidyverse)
library(ggplot2)
library(dplyr)

library(nls2)

# NLME
library(nlme)
library(MplusAutomation)
library(texreg)

# SingleStage
required_packages <- c("mvtnorm", "Matrix", "dplyr", "rjags", 
                       "coda", "tidyr", "lubridate", "tidyverse",
                       "data.table", "nlme")
lapply(required_packages, require, character.only = TRUE)

work_path <- "~/XXYDATAanalysis/GoHiARmodel/substantivedata"
setwd(work_path)

r=2023
C="ECLSK_Fixed_fixedtheta3andMean"
Condi = "Single_Stage_2_large"
dateMark="0709"

set.seed(r)

sigmoid <- function(time, theta1, theta2, theta3){
  theta1*exp(-theta2*exp(-(time)*theta3)) - 3 ### so that the function rage: (-3, +infinite)
}

source("posteriorSummaryStats.R")
 
######################## Data Cleaning ########################
data = read.table("ECLSKpaper.dat")
colnames(data) <- c("CHILDID", "GENDER", "C1_7SC0", "C1R4RSCL", "C1R4RTHT", 
                    "C1R4MSCL", "C1R4MTHT", "C2R4RSCL", "C2R4RTHT", "C2R4MSCL", 
                    "C2R4MTHT", "C3R4RSCL", "C3R4RTHT", "C3R4MSCL", "C3R4MTHT",  
                    "C4R4RSCL", "C4R4RTHT", "C4R4MSCL", "C4R4MTHT",  
                    "C5R4RSCL", "C5R4RTHT", "C5R4MSCL", "C5R4MTHT", 
                    "C6R4RSCL", "C6R4RTHT", "C6R4MSCL", "C6R4MTHT",  
                    "C7R4RSCL", "C7R4RTHT", "C7R4MSCL", "C7R4MTHT",  
                    "WKPARED", "WKSESL", "WKPOV_R", "W1PARED", "W1SESL", 
                    "W3PARED", "W3SESL", "W5PARED", "W5SESL",   
                    "P7NUMSIB", "P7HTOTAL", "W8PARED", "W8SESL", 
                    "C7SPORTS", "C7HRSCLB", "C7HRSRD", "C7TVWKDY", 
                    "C7TVWKEN", "C7VIDWKD", "C7VIDWKN", "C7INTWKD", "C7INTWKN", 
                    "C1read_R", "C1math_R", "C1genT_R", "C2read_R", 
                    "C2math_R", "C2gen_R", "C3read_R", "C3math_R", "C3gen_R", 
                    "C4read_R", "C4math_R", "C4gen_R", "C5read_R", "C5math_R", 
                    "C5gen_R", "C6read_R", "C6math_R", "C6gen_R", "C7read_R", 
                    "C7math_R", "C7gen_R")

data[data==-999] <- NA

data <- mutate(data, 
               read0 = (C1read_R+C2read_R)/2,
               read1 = (C3read_R+C4read_R)/2,
               read3 = C5read_R,
               read5 = C6read_R,
               read8 = C7read_R,
               math0 = (C1math_R+C2math_R)/2,
               math1 = (C3math_R+C4math_R)/2,
               math3 =  C5math_R,
               math5 =  C6math_R,
               math8 =  C7math_R)

ids = unique(data$CHILDID)
time = c(0, 1, 3, 5, 8)
nP=length(ids)
nT=length(time)

data_fit <- matrix(NA,nrow=length(ids)*length(time),ncol=3)
colnames(data_fit) = c("id","time","read")
for( i in 1:nP){
  for(t in 1:nT){
    data_fit[(i-1)*nT+t,"id"]=ids[i]
    data_fit[(i-1)*nT+t,"time"]=time[t]
    data_fit[(i-1)*nT+t,"read"]=data[data$CHILDID==ids[i],paste0("read",time[t])]
  }
}
data_fit = data.frame(data_fit)

#------------------- Missing data analysis -------------------
naniar::vis_miss(data_fit)
naniar::gg_miss_fct(x = data_fit, fct = time)

#------------------- Data visualization -------------------
hist(data_fit$read)

sampleid = sample(ids, 100,replace = FALSE)
(p_1 <- ggplot(data_fit[data_fit$id %in% sampleid,], aes(x = time, group=id)) +
  geom_line(aes(y = read, color="read"),size=0.5) +
  scale_color_manual(values = c("read" = "#8a317e"))+
  labs(title = "read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")+
  scale_x_continuous(breaks = unique(data_fit$time)))

#------------------- Handling missing data -------------------
# remove cases with nT < 5
summary(data_fit$read)

data_fit_noNA = data_fit[complete.cases(data_fit), ]

data_fit_noNA_filtered = data_fit_noNA
for(i in unique(data_fit_noNA_filtered$id)){ #i=1
  dat_i = data_fit_noNA_filtered[which(data_fit_noNA_filtered$id==i),]
  if(nrow(dat_i)<5){data_fit_noNA_filtered=filter(data_fit_noNA_filtered,id!=i)}
}

data_fit_noNA_filtered=as.data.frame(data_fit_noNA_filtered)
sampleid = sample(unique(data_fit_noNA_filtered$id), 100,replace = FALSE)
(p_2 <- ggplot(data_fit[data_fit$id %in% sampleid,], aes(x = time, group=id)) +
    geom_line(aes(y = read, color="read"),size=0.5) +
    scale_color_manual(values = c("read" = "#8a31ee"))+
    labs(title = "read") +
    theme_minimal()+
    theme(axis.title.x = element_text(size = 20), 
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")+
    scale_x_continuous(breaks = unique(data_fit$time)))

# Intercept missing cases so that time are equally spaced
data_fit_equtime <- matrix(NA,
                           nrow=length(unique(data_fit_noNA_filtered$id))*9,
                           ncol=3
                           )
row=1
for(id_now in unique(data_fit_noNA_filtered$id)){
  data_fit_equtime[row:(row+8),1] = id_now #id
  data_fit_equtime[row:(row+8),2] = 0:8 #time
  data_fit_equtime[row+0,3] = filter(data_fit_noNA_filtered, id==id_now, time==0)[,3] #y_t=0
  data_fit_equtime[row+1,3] = filter(data_fit_noNA_filtered, id==id_now, time==1)[,3] #y_t=1
  data_fit_equtime[row+3,3] = filter(data_fit_noNA_filtered, id==id_now, time==3)[,3] #y_t=3
  data_fit_equtime[row+5,3] = filter(data_fit_noNA_filtered, id==id_now, time==5)[,3] #y_t=5
  data_fit_equtime[row+8,3] = filter(data_fit_noNA_filtered, id==id_now, time==8)[,3] #y_t=8
  row=row+9
}
data_fit_equtime = as.data.frame(data_fit_equtime)
colnames(data_fit_equtime)=colnames(data_fit)

sampleid = sample(unique(data_fit_equtime$id), 100,replace = FALSE)
(p_3 <- ggplot(data_fit[data_fit$id %in% sampleid,], aes(x = time, group=id)) +
    geom_line(aes(y = read, color="read"),size=0.5) +
    geom_point(aes(y = read, color="read"),size=1) +
    scale_color_manual(values = c("read" = "orange"))+
    labs(title = "read") +
    theme_minimal()+
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.position = "none")+
    scale_x_continuous(breaks = 0:8))

######################## Descriptives ########################
psych::describe(data_fit_equtime$read)

######################## Starting Points ########################
theta10_start = 4.40
theta20_start = 0.76
theta30_start = 0.58

# Check whether the starting point is reasonable
baseline_data <- matrix(NA,nrow=9,ncol=3)
colnames(baseline_data) = c("time","read","id")
for(t in 0:8){
  baseline_data[t+1,1] = t
  baseline_data[t+1,2] = sigmoid(t,theta10_start,theta20_start,theta30_start)
  baseline_data[t+1,3] = 0
}
baseline_data=as.data.frame(baseline_data)
(p_4 <- ggplot(data_fit_noNA_filtered,
               aes(x = time, group=id)) +
    geom_line(aes(y = read, color="read"),size=0.5) +
    geom_point(aes(y = read, color="read"),size=1) +
    scale_color_manual(values = c("read" = "orange"))+
    labs(title = paste0("start: ",round(theta10_start,2)," | ",
                        round(theta20_start,2)," | ",
                        round(theta30_start,2))) +
    theme_minimal()+
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.position = "none")+
    scale_x_continuous(breaks = 0:8) +
    geom_line(data = baseline_data, aes(x=time, y=read, group=id,colour="#000099"),
              show_guide = FALSE,size=1))

  
baseline_data <- matrix(NA,nrow=9,ncol=3)
colnames(baseline_data) = c("time","read","id")
for(t in 0:8){
  baseline_data[t+1,1] = t
  baseline_data[t+1,2] = sigmoid(t,theta10_start,theta20_start,theta30_start)
  baseline_data[t+1,3] = 0
}
baseline_data=as.data.frame(baseline_data)

######################## Single-Stage Method ########################
file_model = paste0(Condi,"_model","_r",r,"_",dateMark,".txt")
file_ROPEsamples = paste0(Condi,"_ROPEsamples","_r",r,"_",dateMark,".csv")
file_level2pars = paste0(Condi,"_level2pars","_r",r,"_",dateMark,".csv")
file_jagsoutput = paste0(Condi,"_jagsoutput","_r",r,"_",dateMark,".csv")

#------------------- Pre-process data -------------------
outlier_ids <- c("312", "1871")
data_fit_equtime = filter(data_fit_equtime,!(id %in% outlier_ids))

# Get dimensions
ids = unique(data_fit_equtime$id)
ydim=1
nT=9

# Mutate the structure of dat to Y[nT,dim(y+time),id]
Obs_empty = rep(NA, length(ids)*nT*(ydim+1))
Obs = array(Obs_empty,dim=c(nT, (ydim+1), length(ids)))
p=1
for (i in ids){ #i=2
  subset_i = subset(data_fit_equtime,id==i)
  Obs[,,p] = as.matrix(subset_i[,-1])
  p=p+1
} #@@@ it should be noted that Obs[,,p] doesn't correspond to id=p
Obs[,,2]  # data for id=1; columns = (time, y1)

#------------------- MCMC -------------------
# Create a list of all the variables 
jagsData <- list("ids" = ids, "nT" = nT, "Obs" = Obs,
                 "theta10_start" = theta10_start,
                 "theta20_start" = theta20_start,
                 "theta30_start" = theta30_start)

# Specify the Structured Latent Curve Model
SLCmodel_0 = cat("
model {
# PART 1. Specifying Level-1 model
  for (id in 1:length(ids)){ # opening loop for person (id)

    Obs[1,2,id]~dnorm(MU[id,1],1/IIV)
    MU[id,1]<-MU_thetas[id,1]*exp(-MU_thetas[id,2]*exp(-Obs[1,1,id]*MU_thetas[id,3]))-3
    
    for (t in 2:nT) { 
      Obs[t,2,id] ~ dnorm(MU[id,t] + AR[id]*(Obs[t-1,2,id]-MU[id,t-1]), 1/IIV)
      MU[id,t] <- MU_thetas[id,1]*exp(-MU_thetas[id,2]*exp(-Obs[t,1,id]*MU_thetas[id,3]))-3
    } 
    
# PART 2. Specifying Level-2 model # should match data generation
  AR[id] ~ dnorm(Level2MeanAR, 1/Level2Var_AR)T(-1,1)
  MU_thetas[id,1] ~ dnorm(Level2Mean_theta1, 1/Level2Var_theta1)
  MU_thetas[id,2] ~ dnorm(Level2Mean_theta2, 1/Level2Var_theta2)
  MU_thetas[id,3] ~ dnorm(Level2Mean_theta3, 1/Level2Var_theta3)
    
 } # closing loop for  person (id)       
      
# PART 3. Specifying prior distributions 
  # AR
    Level2MeanAR ~ dunif(-0.5,0.5)
    Level2Sd_AR ~ dunif(0,0.5)
    Level2Var_AR = pow(Level2Sd_AR,2)
  # MU
    Level2Mean_theta1 ~ dunif(3.5,5.5)
    Level2Mean_theta2 ~ dunif(0,1.5)
    Level2Mean_theta3 = 0.596
    
    Level2Sd_theta1 ~ dunif(0,0.9)
    Level2Sd_theta2 ~ dunif(0,0.4)
    #log_Level2Sd_theta3 ~ dunif(-6.9077553,-0.9162907) #log(c(0.001,0.4))
    #Level2Sd_theta3 = exp(log_Level2Sd_theta3)
    #Level2Sd_theta3 ~ dunif(0,0.4)
    
    Level2Var_theta1 = pow(Level2Sd_theta1,2)
    Level2Var_theta2 = pow(Level2Sd_theta2,2)
    Level2Var_theta3 = 0.0001
    
  # IIV # IIV = 1
    IIV <- exp(logIIV)   
    logIIV ~ dunif(-13.81551, 0) #log(c(0.001^2,1^2))
    #IIV ~ dunif(0,0.5)
}
",file = file_model)

# parameters to track
summ_rownames = c("IIV",
                  "Level2MeanAR","Level2Var_AR",
                  "Level2Mean_theta1","Level2Mean_theta2","Level2Mean_theta3",
                  "Level2Var_theta1","Level2Var_theta2","Level2Var_theta3")
parameters <- c(summ_rownames)#"AR","MU","MU_thetas",

# Specifying sampler settings
adaptation  <- 5000
burnin  <- 5000
chains  <- 2
thinning  <- 5
# Defining the number of posterior samples per chain for JAGS
nrOfIter <- 50000

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

######################## Results Diagnosis ########################
# Check trace plots
#codaSamples_2 <- lapply(codaSamples, function(x) x[, 2144:2152])
#str(codaSamples_2)
#coda::traceplot(codaSamples_2)
coda::traceplot(codaSamples)
