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
C="ECLSK"
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

######################## Searching Starting Points ########################
theta10_start = 4.44
theta20_start = 0.67
theta30_start = 0.42

# theta10_start = 4.40
# theta20_start = 0.76
# theta30_start = 0.58

baseline_data <- matrix(NA,nrow=9,ncol=3)
colnames(baseline_data) = c("time","read","id")
for(t in 0:8){
  baseline_data[t+1,1] = t
  baseline_data[t+1,2] = sigmoid(t,theta10_start,theta20_start,theta30_start)
  baseline_data[t+1,3] = 0
}
baseline_data=as.data.frame(baseline_data)

# Check whether the starting point is reasonable
(p_4 <- ggplot(data_fit_noNA_filtered[data_fit_noNA_filtered$id %in% sampleid,],
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

#------------------- NLSLIST -------------------
# Fit a Gompertz curve for each individual to find the range of parameters
data_fit_NLSLIST = data_fit_noNA_filtered

# Get dimensions
ids = unique(data_fit_NLSLIST$id)

# Detrend Data
read_sigmoid.nlsList = nlsList(read~sigmoid(time, theta1, theta2, theta3)|id,
                               data=data_fit_NLSLIST,
                               start = c(theta1 = theta10_start, # NLME
                                         theta2 = theta20_start,
                                         theta3 = theta30_start),
                               control = nls.control(maxiter = 10000, 
                                                     minFactor = 1e-10,
                                                     warnOnly = T, 
                                                     tol=1e-4)
)

thetas_pre = coefficients(read_sigmoid.nlsList)

dat_detr <- matrix(NA, nrow=nrow(data_fit_NLSLIST), ncol=ncol(data_fit_NLSLIST))
colnames(dat_detr) <- colnames(data_fit_NLSLIST)

dat_detr[,c("id","time")] = as.matrix(data_fit_NLSLIST[,c("id","time")])

thetas_pre = coefficients(read_sigmoid.nlsList)

unitime=unique(data_fit_NLSLIST$time)

for(i in 1:length(unique(data_fit_NLSLIST$id))){ # i=1
  id_now=unique(data_fit_NLSLIST$id)[i]
  theta_i = thetas_pre[i,]
  y_true_i = filter(data_fit_NLSLIST, id==id_now)$read
  if(sum(is.na(theta_i))>0){
    # dat_i = y_true_i
    # dat_detr[((i-1)*nT+1):(i*nT),"read"]=as.matrix(dat_i)
    dat_detr[((i-1)*nT+1):(i*nT),"read"]=rep(NA,5)
  }else{
    for(t in 1:nT){
      y_hat_i = sigmoid(unitime[t], theta_i[1], theta_i[2], theta_i[3])
      dat_i = y_true_i[t] - y_hat_i
      if(abs(dat_i)>5){
        dat_detr[(i-1)*nT+t,"read"] <- NA
      }else{dat_detr[(i-1)*nT+t,"read"]=as.matrix(dat_i)}
    }
  }
}

dat_detr <- as.data.frame(dat_detr)

for(i in unique(dat_detr$id)){ #i=1
  dat_i = dat_detr[which(dat_detr$id==i),]
  if(sum(is.na(dat_i$read))>0){dat_detr=filter(dat_detr,id!=i)}
}

length(unique(dat_detr$id))/length(unique(data_fit_NLSLIST$id))

# Check Data
IDtoPlot = sort(sample(unique(dat_detr$id), 20, replace = FALSE))
data_subset_1 = data_fit_NLSLIST[data_fit_NLSLIST$id %in% IDtoPlot,]
data_subset_1$group <- "trended"
data_subset_1_detrended =dat_detr[dat_detr$id %in% IDtoPlot,]
data_subset_1_detrended$group="detrended"

y_limits <- range(data_subset_1$read)
x_limits <- range(data_subset_1$time)
#x_breaks <- seq(x_limits[1], x_limits[2], by = 1)


#------------------- Check parameters & Remove outliers -------------------
psych::describe(thetas_pre)
boxplot(thetas_pre$theta1)
boxplot(thetas_pre$theta2)
boxplot(thetas_pre$theta3)

plot(thetas_pre$theta1)
plot(thetas_pre$theta2)
plot(thetas_pre$theta3)

# According to plots, there exist two outliers
# Let's check these two samples row=312, row=1871
# rownames(thetas_pre) <- NULL
outlier_ids <- c("312", "1871")
thetas_pre[outlier_ids,]

data_subset_outlier <- filter(data_fit_NLSLIST, id %in% outlier_ids)
ggplot(data_subset_outlier, aes(x = time, y = read, group = id)) +
  geom_point(size=0.5, color ="#00BFC4") +
  geom_line(color ="#00BFC4") +
  scale_x_continuous(breaks = unique(data_subset_1$time))+
  scale_y_continuous(limits = y_limits) +
  ggtitle("Outliers: trended individual traces") +
  xlab("Time") +
  ylab("read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_blank(),
        legend.position = "none")

ggplot(data_subset_outlier, aes(x = time, y = read, group = id)) +
  geom_point(size=0.5, color="#F8766D") +
  geom_line(color="#F8766D") +
  scale_x_continuous(breaks = unique(data_subset_1$time))+
  scale_y_continuous(limits = y_limits) +
  ggtitle("Outliers:detrended individual traces") +
  xlab("Time") +
  ylab("read") +
  theme_minimal()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_blank(),
        legend.position = "none")

# Remove these two students
data_fit_NLSLIST_rmoutliers <- filter(data_fit_NLSLIST,!(id %in% outlier_ids))

length(unique(data_fit_NLSLIST_rmoutliers$id))

#------------------- Fit NLSLIST again -------------------
# Fit a Gompertz curve for each individual to find the range of parameters
data_fit_NLSLIST = data_fit_NLSLIST_rmoutliers

# Get dimensions
ids = unique(data_fit_NLSLIST$id)

# Fit Data
read_sigmoid.nlsList = nlsList(read~sigmoid(time, theta1, theta2, theta3)|id,
                               data=data_fit_NLSLIST,
                               start = c(theta1 = theta10_start, # NLME
                                         theta2 = theta20_start,
                                         theta3 = theta30_start),
                               control = nls.control(maxiter = 10000, 
                                                     minFactor = 1e-10,
                                                     warnOnly = T, 
                                                     tol=1e-4)
)

thetas_pre = coefficients(read_sigmoid.nlsList)
res_detr = summary(read_sigmoid.nlsList)
# Check model fitting 
residuals <- residuals(read_sigmoid.nlsList)
plot(residuals)

#------------------- Check parameters -------------------
psych::describe(thetas_pre)

boxplot(thetas_pre$theta1)
boxplot(thetas_pre$theta2)
boxplot(thetas_pre$theta3)

plot(thetas_pre$theta1)
plot(thetas_pre$theta2)
plot(thetas_pre$theta3)

#------------------- Detrend parameters -------------------
dat_detr <- matrix(NA, nrow=nrow(data_fit_NLSLIST), ncol=ncol(data_fit_NLSLIST))
colnames(dat_detr) <- colnames(data_fit_NLSLIST)

dat_detr[,c("id","time")] = as.matrix(data_fit_NLSLIST[,c("id","time")])

unitime=unique(data_fit_NLSLIST$time)

for(i in 1:length(unique(data_fit_NLSLIST$id))){ # i=1
  id_now=unique(data_fit_NLSLIST$id)[i]
  theta_i = thetas_pre[i,]
  y_true_i = filter(data_fit_NLSLIST, id==id_now)$read
  if(sum(is.na(theta_i))>0){
    dat_detr[((i-1)*nT+1):(i*nT),"read"]=rep(NA,5)
  }else{
    for(t in 1:nT){
      y_hat_i = sigmoid(unitime[t], theta_i[1], theta_i[2], theta_i[3])
      dat_i = y_true_i[t] - y_hat_i
      if(abs(dat_i)>5){
        dat_detr[(i-1)*nT+t,"read"] <- NA
      }else{dat_detr[(i-1)*nT+t,"read"]=as.matrix(dat_i)}
    }
  }
}

dat_detr <- as.data.frame(dat_detr)

for(i in unique(dat_detr$id)){ #i=1
  dat_i = dat_detr[which(dat_detr$id==i),]
  if(sum(is.na(dat_i$read))>0){dat_detr=filter(dat_detr,id!=i)}
}

length(unique(dat_detr$id))/length(unique(data_fit_NLSLIST$id))

# Check Data
IDtoPlot = sort(sample(unique(dat_detr$id), 20, replace = FALSE))
data_subset_1 = data_fit_NLSLIST[data_fit_NLSLIST$id %in% IDtoPlot,]
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
Condi = "NLSLIST"

nP = length(unique(dat_detr_long$id))
  
file_detrdat_long = paste0(Condi,"_detrdat_long_nT",nT,
                           "_nP",nP,"_r",r,"_",dateMark,".csv")
file_model = paste0(Condi,"_model_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".rds")


print(paste0("file_detrdat_long: ", file_detrdat_long))
print(paste0("file_model: ", file_model))


saveRDS(res_detr, file = file_model)
#saveRDS(dat_detr_long, file = file_model)
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
data_diagnosis = data_mplus_equtime
data_mplus_equtime[is.na(data_mplus_equtime)] <- -99999

#------------------- Fit DSEM in Mplus -------------------
# Set Mplus input File
input = mplusObject(TITLE = "A two-level DSEM model for one continuous dependent
                      variable with random intercepts and random slopes",
                    VARIABLE = "MISSING ARE ALL (-99999);
                      CLUSTER = id;
                      LAGGED =  read(1);",
                    ANALYSIS = "TYPE = TWOLEVEL RANDOM;
                      ESTIMATOR = BAYES;
                      ALGORITHM = GIBBS;
                      PROCESSORS = 2;
                      BITERATIONS = 60000(10000);",
                    MODEL = "%WITHIN%
                      s1 | read ON read&1;
                     %BETWEEN%
                    read@0;
                    [read@0];",# not equally spaced 
                    OUTPUT = "TECH1 TECH8;", #STANDARDIZED (CLUSTER)
                    SAVEDATA = paste0("SAVE=FSCORES(50 10);FILE IS ",Condi,"_FSCORES",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".csv;"),
                    PLOT = "TYPE = PLOT3;
                      FACTOR = ALL;",
                    usevariables = c("id","read"),
                    rdata = data_mplus_equtime,
                    autov = TRUE)
res = mplusModeler(input,
                   modelout = paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".inp"), 
                   hashfilename=FALSE,
                   run = 1L)
alloutput=readModels(paste0(Condi,"_",C,"_nT",nT,"_nP",nP,"_r",r,"_",dateMark,".out"))
alloutput$parameters$unstandardized

#------------------- Diagnosis of Residuals -------------------
library(aTSA)
TestTable <- matrix(NA,
                    nrow=length(unique(data_diagnosis$id)),
                    ncol=1+12
)
for (i in 1:length(unique(data_diagnosis$id))){
  # ADF test
  adf_i = adf.test(data_diagnosis[((i-1)*9+1):(i*9),"read"])
  # KPSS test
  kpss_i = kpss.test(data_diagnosis[((i-1)*9+1):(i*9),"read"])
  TestTable[i,1]=data_diagnosis[((i-1)*9+1),"id"]
  
  TestTable[i,2:3]=adf_i$type1[1,2:3]
  TestTable[i,4:5]=adf_i$type2[1,2:3]
  TestTable[i,6:7]=adf_i$type3[1,2:3]
  
  TestTable[i,8:9]=kpss_i[1,2:3]
  TestTable[i,10:11]=kpss_i[2,2:3]
  TestTable[i,12:13]=kpss_i[3,2:3]
}
colnames(TestTable)=c("id","ADF_type1","p_ADF_type1",
                      "ADF_type2","p_ADF_type2",
                      "ADF_type3","p_ADF_type3",
                      "KPSS_type1","p_KPSS_type1",
                      "KPSS_type2","p_KPSS_type2",
                      "KPSS_type3","p_KPSS_type3")

TestTable <- as.data.frame(TestTable) 

# Create an empty data frame to store the results
summary_df <- data.frame()

# Loop through all the types
types <- c("ADF_type1", "ADF_type2", "ADF_type3",
           "KPSS_type1", "KPSS_type2", "KPSS_type3")

for (type in types) {
  # Calculate the mean for each type
  mean_column <- mean(TestTable[, type], na.rm = TRUE)
  summary_df[1,type] = mean_column
  
  # Calculate the proportion for the type
  p_type = paste0("p_",type)
  proportion_column <- sum(TestTable[, p_type] <= 0.01, na.rm = TRUE) / nrow(TestTable)
  summary_df[2,type] <- proportion_column
}
summary_df[3,1:3]=summary_df[2,1:3]
summary_df[3,4:6]=1-summary_df[2,4:6]
rownames(summary_df) <- c("mean", "proportaion_of_p_leq_0_01","proportaion_of_stationary_time_series")

# Display the summarized information
print(summary_df)

# min(summary_df[3,1:3])
# max(summary_df[3,1:3])
# min(summary_df[3,4:6])
# max(summary_df[3,4:6])

round(min(summary_df[3,1:3]),3)
round(max(summary_df[3,1:3]),3)
round(min(summary_df[3,4:6]),3)
round(max(summary_df[3,4:6]),3)


