# This script is for paper 3. We look at the associations between indices of weight instability (WV and WC)
# and associations with cardiometabolic health markers and body composition
# we will use 2 data sets: (1) a WV and WC calculated DF and (2) the main cid data

library(dplyr)
library(lubridate)
library("birk")
library(broom)
library("imputeMulti")
library("imputeTS")
library(zoo)
library(xts)
library(reshape2)
library(plyr)
library("randomForest")
library("Metrics")
library(Hmisc)
library(forecast)
library(Amelia)
library("tidyr")
library('rlist')
library(psych)
library('sjstats')
library(stats)
library('multilevel')
library('gam')
library('ggthemes')
library("forcats")
library('BlandAltmanLeh')
library('anytime')
library('tourr')
library('pracma')
library('rowr')
library(tableone)
library('lme4')

setwd("C:/Users/jaket/Dropbox/PhD/NoHoW Analyses/Weight variability")

rmse_sum<-read.csv('imputation_rmse_summary_table.csv')
rmse_sum[3:4]<-round(rmse_sum[3:4], 2)
rmse_sum$rmse_se<-paste(rmse_sum$error, ' (', rmse_sum$se, ')', sep='')
rmse_sum2<-rmse_sum[,-c(3:4)]
rmse_sum2_wide<-reshape(rmse_sum2, idvar=c('insertion', 'missingness'), direction='w', timevar = 'strategy')
colnames(rmse_sum2_wide)[3:12]<-c('KNN', 'PMM', 'RF', 'Arima Kal', 'EWMA', 
                                  'Linear Interpolation', 'Spline Interpolation', 'Stine Interpolation',
                                  'SMKS', 'TS Clean')
write.csv(rmse_sum2_wide, 'RMSE_results_final_230120.csv', row.names = F)



WV_sum<-read.csv('WV_estimation_summaries.csv')
WV_sum[5:6]<-round(WV_sum[5:6],1)
WV_sum$imputation_strategy<- as.character(WV_sum$imputation_strategy)
WV_sum$imputation_strategy[is.na(WV_sum$imputation_strategy)]<-'None'
WV_sum$mean_se<-paste(WV_sum$mean_deviation, ' (', WV_sum$se, ')', sep='')
WV_sum2<-WV_sum[,-c(2,5:6)]
WV_sum2_long<-reshape(WV_sum2, idvar=c('WV', 'missingness'), direction='w', timevar = 'imputation_strategy')

colnames(WV_sum2_long)[3:13]<-c('None', "Arima Kal","EWMA","KNN","Linear Interpolation",
                                "PMM" ,"RF","Spline Interpolation","Stine Interpolation",
                                "StrucModel Kal" , "TS Clean")
write.csv(WV_sum2_long, 'WV_results_final_230120.csv', row.names = F)

