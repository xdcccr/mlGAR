rm(list=ls())

# Data Cleaning
required_packages <- c("mvtnorm", "Matrix", "dplyr", "tidyverse",
                       "data.table", "forestmangr", "MplusAutomation",
                       "texreg","forestmangr")
lapply(required_packages, require, character.only = TRUE)
# work_path <- "~/XXYDATAanalysis/GoHiARmodel/substantivedata"
work_path <- "/Users/xiaoyuexiong/XXYDATAanalysis/IP1-DetrendingPanelTimeSeries/CodesAndResultsForRevision/RevisedEmpirical"
setwd(work_path)

r=2023
C="ECLSK"
dateMark="2501final"

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

# Setting 
Condi = "Linear"

file_detrdat_long = paste0(Condi,"_detrdat_long_nT",nT,
                           "_nP",nP,"_r",r,"_",dateMark,".csv")
file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")


print(paste0("file_detrdat_long: ", file_detrdat_long))
print(paste0("file_model: ", file_model))

######################### Two-Stage Linear ########################
data_fit_Linear = data_fit_noNA_filtered

outlier_ids <- c("312", "1871")

data_subset_outlier <- filter(data_fit_Linear, id %in% outlier_ids)
ggplot(data_subset_outlier, aes(x = time, y = read, group = id)) +
  geom_point(size=0.5, color ="#00BFC4") +
  geom_line(color ="#00BFC4") +
  scale_x_continuous(breaks = unique(data_fit_Linear$time))+
  #scale_y_continuous(limits = y_limits) +
  ggtitle("Outliers: trended individual traces") +
  xlab("Time") +
  ylab("read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_blank(),
        legend.position = "none")

data_fit_Linear = filter(data_fit_noNA_filtered, !(id %in% outlier_ids))
#------------------- Fit Data -------------------
# Get dimensions
ids = unique(data_fit_Linear$id)

# Detrend Data
read_Linear = forestmangr::lm_table(data_fit_Linear, read ~ time, "id")

# 计算N1-Linear的总统计量
n_subjects <- 2143  # 个体数量
n_timepoints <- 5   # 每个个体的时间点数量
n_total <- n_subjects * n_timepoints  # 总样本量 = 10,715

k_per_model <- 3  # 每个线性模型的参数数量（截距、斜率、残差方差）
k_total <- k_per_model * n_subjects  # 总参数数量 = 6,429

# 计算总deviance
total_deviance <- 0
for (i in unique(data_fit_Linear$id)) {
  id_data <- subset(data_fit_Linear, id == i)
  model <- lm(read ~ time, data = id_data)
  total_deviance <- total_deviance + deviance(model)
}

# 基于总deviance计算BIC
total_BIC <- total_deviance + k_total * log(n_total)

# 打印结果
cat("Total deviance:", total_deviance, "\n")
cat("Total parameters:", k_total, "\n")
cat("Total observations:", n_total, "\n")
cat("Total BIC:", total_BIC, "\n")

dat_detr <- matrix(NA, nrow=nrow(data_fit_Linear), 
                   ncol=ncol(data_fit_Linear))

colnames(dat_detr) <- colnames(data_fit_Linear)

dat_detr[,c("id","time")]=as.matrix(data_fit_Linear[,c("id","time")])

for(i in 1:length(unique(data_fit_Linear$id))){
  b0_i=read_Linear$b0[i]
  b1_i=read_Linear$b1[i]
  for(t in 1:nT){
    dat_detr[(i-1)*nT+t,"read"] = data_fit_Linear[(i-1)*nT+t,"read"]-b0_i-b1_i*dat_detr[(i-1)*nT+t,"time"] 
  }
}

dat_detr <- as.data.frame(dat_detr)

# Check Data
IDtoPlot = sort(sample(1:nP, 20, replace = FALSE))
data_subset_1 = data_fit_Linear[data_fit_Linear$id %in% IDtoPlot,]
data_subset_1$group <- "trended"
data_subset_1_detrended =dat_detr[dat_detr$id %in% IDtoPlot,]
data_subset_1_detrended$group="detrended"

y_limits <- range(data_subset_1$read)
x_limits <- range(data_subset_1$time)
#x_breaks <- seq(x_limits[1], x_limits[2], by = 1)

ggplot(data_subset_1, aes(x = time, y = read, group = id)) +
  geom_point(size=0.5, color ="#00BFC4") +
  geom_line(color ="#00BFC4") +
  scale_x_continuous(breaks = unique(data_subset_1$time))+
  scale_y_continuous(limits = y_limits) +
  ggtitle("trended individual traces") +
  xlab("Time") +
  ylab("read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_blank(),
        legend.position = "none")

ggplot(data_subset_1_detrended, aes(x = time, y = read, group = id)) +
  geom_point(size=0.5, color="#F8766D") +
  geom_line(color="#F8766D") +
  scale_x_continuous(breaks = unique(data_subset_1$time))+
  scale_y_continuous(limits = y_limits) +
  ggtitle("detrended individual traces") +
  xlab("Time") +
  ylab("read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_blank(),
        legend.position = "none")

dat_detr_long <- as.data.frame(dat_detr)

# Save data
saveRDS(dat_detr, file = file_model)
write.table(dat_detr_long, file_detrdat_long, row.names = FALSE, col.names = TRUE, quote = F)

dat_mplus = dat_detr_long

# Intercept missing cases so that time are equally spaced
data_mplus_equtime <- matrix(NA,
                             nrow=length(unique(dat_mplus$id))*9,
                             ncol=3
)
row=1
for(id_now in unique(dat_mplus$id)){
  data_mplus_equtime[row:(row+8),1] = id_now #id
  data_mplus_equtime[row:(row+8),2] = 0:8 #time
  data_mplus_equtime[row+0,3] = filter(dat_mplus, id==id_now, time==0)[,3] #y_t=0
  data_mplus_equtime[row+1,3] = filter(dat_mplus, id==id_now, time==1)[,3] #y_t=1
  data_mplus_equtime[row+3,3] = filter(dat_mplus, id==id_now, time==3)[,3] #y_t=3
  data_mplus_equtime[row+5,3] = filter(dat_mplus, id==id_now, time==5)[,3] #y_t=5
  data_mplus_equtime[row+8,3] = filter(dat_mplus, id==id_now, time==8)[,3] #y_t=8
  row=row+9
}
data_mplus_equtime = as.data.frame(data_mplus_equtime)
colnames(data_mplus_equtime)=colnames(dat_mplus)
# data_diagnosis = data_mplus_equtime
# data_mplus_equtime[is.na(data_mplus_equtime)] <- -99999

#------------------- Fit DSEM in rJAGS -------------------
data_JAGS_equtime =  data_mplus_equtime
file_model = paste0(C,"_",Condi,"_model_nT",nT,"_nP",nP,"_",dateMark,".txt")
file_ROPEsamples = paste0(C,"_",Condi,"_ROPEsamples_nT",nT,"_nP",nP,"_",dateMark,".csv")
file_level2pars = paste0(C,"_",Condi,"_level2pars_nT",nT,"_nP",nP,"_",dateMark,".csv")
file_jagsoutput = paste0(C,"_",Condi,"_jagsoutput_nT",nT,"_nP",nP,"_",dateMark,".csv")

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
level2_params <- resulttable[c("IIV","Level2MeanAR","Level2Var_AR","Deviance"),]
write.csv(level2_params, file=file_level2pars)

# Save full JAGS output
sumjags <- summary(codaSamples)
write.csv(sumjags$statistics, file=file_jagsoutput)
