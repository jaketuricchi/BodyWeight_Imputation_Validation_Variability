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
library('RColorBrewer')
library(tableone)
library(anytime)
library(simsem)


#load df
#setwd("~/Dropbox/PhD/NoHoW Analyses/Weight variability")
setwd("C:/Users/jaket/Dropbox/PhD/NoHoW Analyses/Weight variability")

#load df with true weights
df_true<-read.csv('most_complete5.csv')%>%arrange(ID, day_no)%>%subset(select=c(weight, ID, day_no))

#load df with simulation missingness
df_rpm<-read.csv('rpm_df.csv')%>%subset(select=c(-weight))
df_mcar<-read.csv('NA_df_50_MCAR_161219.csv')
df_mcar$day_no=df_true$day_no
df_mcar<-df_mcar[,c(2:82,1)]

#load imputed dfs
df_imp_mcar_uni<-read.csv('basic_imp_big.csv')%>%subset(select=c(-X))
df_imp_mcar_uni$day_no<-df_true$day_no
df_imp_mcar_uni<-df_imp_mcar_uni[,c(2:561, 1, 562)]
df_imp_rpm_uni<-read.csv('RPM_univariate_imp_df.csv')
df_imp_rpm_uni$ID<-df_true$ID
df_imp_rpm_uni$day_no<-df_true$day_no
df_imp_mcar_multi<-read.csv('multi_imp_mcar_df.csv')
df_imp_mcar_multi$ID<-df_true$ID
df_imp_mcar_multi$day_no<-df_true$day_no
df_imp_rpm_multi<-read.csv('multi_imp_rpm_df.csv')
df_imp_rpm_multi$ID<-df_true$ID
df_imp_rpm_multi$day_no<-df_true$day_no

#merge all
dfs<-list(df_true, df_mcar, df_rpm, df_imp_mcar_multi, df_imp_mcar_uni, df_imp_rpm_multi, df_imp_rpm_uni)

#calculate linear RMSE for all dfs
rmse_rel_fn<-function(df){
  colnames(df)[1]<-'ts'
  res<-residuals(lm (ts~day_no, data=df))
  weights<-na.omit(df$ts)
  res<-cbind.data.frame(res, weights)
  rel_res<-(res$res/res$weights)*100
  RMSE=data.frame(sqrt(mean(rel_res^2)))
  colnames(RMSE)[1]<-'RMSE'
  return(RMSE)
}

#loop over each col
out_list<-list()
multi_loop_rmse_fn<-function(x){
  columns_n<-ncol(x)-2
  for (i in 1:columns_n){
    df<-data.frame(x[, c(i)])
    df$day_no<-x$day_no
    rmse<-rmse_rel_fn(df)
    colnames(rmse)<-colnames(x)[i]
    out_list[[i]]<-rmse
  }
  out<-bind_cols(out_list)
  return(out)
}

#split function by ID
split_apply_rmse_fn<-function(x){
  rmse_list<-dlply(x, 'ID', multi_loop_rmse_fn)
  rmse_df<-bind_rows(rmse_list)
  return(rmse_df)
}

#do this for all dfs in list
rmse_lists<-lapply(dfs, split_apply_rmse_fn)

#add IDs to each list
IDs<-data.frame(levels(df_true$ID))
colnames(IDs)<-'ID'

rmse_lists2<-lapply(rmse_lists, function(x) 
  cbind(x, ID = IDs, WV='RMSE'))

#melt all
rmse_lists_melt<-lapply(rmse_lists2, function(x) melt(x, id.var=c('ID', 'WV')))


##############################
## NLMD #####################
##############################
######################################proof of concept working##################################################
x2<-dfs[[1]]
x<-filter(x2, ID=='CPH_HAP_8308')
x<-filter(x, !is.na(weight))
weights<-as.matrix(x$weight)
x$na10<-imposeMissing(data.mat=weights, pmMCAR=0.1)
x$na20<-imposeMissing(data.mat=weights, pmMCAR=0.2)
x$na30<-imposeMissing(data.mat=weights, pmMCAR=0.3)
x$na40<-imposeMissing(data.mat=weights, pmMCAR=0.4)
x$na50<-imposeMissing(data.mat=weights, pmMCAR=0.5)
x$na60<-imposeMissing(data.mat=weights, pmMCAR=0.6)
x$na70<-imposeMissing(data.mat=weights, pmMCAR=0.7)
x$na80<-imposeMissing(data.mat=weights, pmMCAR=0.8)
x$na90<-imposeMissing(data.mat=weights, pmMCAR=0.9)
x$loess_weights5<-loess(x$na50~x$day_no,span=0.985-0.0014*n_weights)$fit #correct
ggplot(x, aes(x=day_no, y=weight))+geom_point()+geom_point(aes(y=loess_weights5))
n_weights<-as.numeric(sum(!is.na(x$weight)))
NLMD_correct<-mean(abs(loess(x$weight~x$day_no,span=(0.985-0.0014*n_weights)$residuals)))
NLMD_10<-mean(abs(loess(x$na10~x$day_no,span=0.6)$residuals))
NLMD_20<-mean(abs(loess(x$na20~x$day_no,span=0.61)$residuals)) 
NLMD_30<-mean(abs(loess(x$na30~x$day_no,span=0.62)$residuals)) 
NLMD_40<-mean(abs(loess(x$na40~x$day_no,span=0.65)$residuals)) 
NLMD_50<-mean(abs(loess(x$na50~x$day_no,span=0.67)$residuals)) 
NLMD_60<-mean(abs(loess(x$na60~x$day_no,span=0.7)$residuals)) 
NLMD_70<-mean(abs(loess(x$na70~x$day_no,span=0.85)$residuals)) 
NLMD_80<-mean(abs(loess(x$na80~x$day_no,span=0.95)$residuals)) 
NLMD_90<-mean(abs(loess(x$na90~x$day_no,span=1)$residuals)) 
xx<-c(sum(!is.na(x$weight)), sum(!is.na(x$na10)), sum(!is.na(x$na20)), sum(!is.na(x$na30)), sum(!is.na(x$na40)),
      sum(!is.na(x$na50)), sum(!is.na(x$na60)), sum(!is.na(x$na70)), sum(!is.na(x$na80)), sum(!is.na(x$na90)))
y<-c(0.55, 0.6, 0.61, 0.62, 0.65, 0.67, 0.7, 0.85, 0.95, 1)
trend<-as.data.frame(cbind(xx, y))
summary(lm<-lm(y~xx))
coef(lm)
ggplot(data=trend, aes(x=xx, y=y))+geom_point()+geom_smooth(method='lm', se=F)
########################################################################################################################

NLMD_fn<-function(df){
  colnames(df)[1]<-'ts'
  loess_data<-subset(df, select=c(ts, day_no))
  loess_data=na.omit(loess_data)
  n_weights<-sum(!is.na(df$weight))
  loess_residuals<-loess(loess_data$ts~loess_data$day_no,span=0.985-0.0014*n_weights)$residuals
  loess_rel_res<-(loess_residuals/loess_data$ts)*100
  NLMD<-data.frame(mean(abs(loess_rel_res)))
  colnames(NLMD)<-'NLMD'
  return(NLMD)
}


#loop over each col
out_list_nlmd<-list()
multi_loop_NLMD_fn<-function(x){
  columns_n<-ncol(x)-2
  for (i in 1:columns_n){
    df<-data.frame(x[, c(i)])
    df$day_no<-x$day_no
    NLMD<-NLMD_fn(df)
    colnames(NLMD)<-colnames(x)[i]
    out_list_nlmd[[i]]<-NLMD
  }
  out<-bind_cols(out_list_nlmd)
  return(out)
}

#split function by ID
split_apply_nlmd_fn<-function(x){
  nlmd_list<-dlply(x, 'ID', multi_loop_NLMD_fn)
  nlmd_df<-bind_rows(nlmd_list)
  return(nlmd_df)
}

#do this for all dfs in list
nlmd_lists<-lapply(dfs, split_apply_nlmd_fn)

#add IDs
nlmd_lists2<-lapply(nlmd_lists, function(x) 
  cbind(x, ID = IDs, WV='NLMD'))

#melt all
nlmd_lists_melt<-lapply(nlmd_lists2, function(x) melt(x, id.var=c('ID', 'WV')))

#extract information from colnames on missingness, imputation etc
#missingness

extract_info_fn<-function(x){
  #how much simulated missing data?
  x$missingness<-as.factor(ifelse(grepl('_20_', x$variable), '20',
                           ifelse(grepl('X20_', x$variable), '20',
                           ifelse(grepl('_40_', x$variable), '40',
                                  ifelse(grepl('X40_', x$variable), '40',
                                  ifelse(grepl('_60_', x$variable), '60',
                                         ifelse(grepl('X60_', x$variable), '60',
                                         ifelse(grepl('_80_', x$variable), '80',
                                                ifelse(grepl('X80_', x$variable), '80', '0')))))))))
  
  #missing/imputed
  x$missing_imputed<-as.factor(ifelse(grepl('_20_', x$variable), 'imputed',
                                      ifelse(grepl('X20_', x$variable), 'missing',
                                             ifelse(grepl('_40_', x$variable), 'imputed',
                                                    ifelse(grepl('X40_', x$variable), 'missing',
                                                           ifelse(grepl('_60_', x$variable), 'imputed',
                                                                  ifelse(grepl('X60_', x$variable), 'missing',
                                                                         ifelse(grepl('_80_', x$variable), 'imputed',
                                                                                ifelse(grepl('X80_', x$variable), 'missing', 'real')))))))))
  
  
  #imputation
  x$imputation_strategy<-ifelse(grepl("TsC_", x$variable), 'TS Clean',
                                ifelse(grepl("Lint_", x$variable), 'Linear Interpolation',
                                       ifelse(grepl("Stine_", x$variable), 'Stine Interpolation',
                                              ifelse(grepl("Spline_", x$variable), 'Spline Interpolation',
                                                     ifelse(grepl("StrKal_", x$variable), 'StrucModel Kal',
                                                            ifelse(grepl("AutoKal_", x$variable), 'Arima Kal',
                                                                   ifelse(grepl("EWMA_", x$variable), 'Exp Weighted MA',
                                                                          ifelse(grepl("KNN",x$variable), 'KNN', 
                                                                            ifelse(grepl("PMM",x$variable), 'PMM',
                                                                                ifelse(grepl("RF",x$variable), 'RF', NA))))))))))
  
  
  #simulation
  x$simulation<-sub(".*_", "", x$variable)
  
  return(x)
}

rmse_lists3<-lapply(rmse_lists_melt, extract_info_fn)
nlmd_lists3<-lapply(nlmd_lists_melt, extract_info_fn)

#extract true nlmd and rmse
truth_rmse<-rmse_lists3[[1]]%>%subset(select=c(ID, WV, value))
colnames(truth_rmse)[3]<-'true_value'
truth_nlmd<-nlmd_lists3[[1]]%>%subset(select=c(ID, WV, value))
colnames(truth_nlmd)[3]<-'true_value'

calculate_WV_diff_rmse<-function(x){
  both_wv<-merge(x, truth_rmse, by=c('ID','WV'), all=T)%>%
    mutate(relative_diff=100-(true_value/value)*100)
  return(both_wv)
}

calculate_WV_diff_nlmd<-function(x){
  both_wv<-merge(x, truth_nlmd, by=c('ID','WV'), all=T)%>%
    mutate(relative_diff=100-(true_value/value)*100)
  return(both_wv)
}

rmse_list_with_diff<-lapply(rmse_lists3[2:7], calculate_WV_diff_rmse)
nlmd_list_with_diff<-lapply(nlmd_lists3[2:7], calculate_WV_diff_nlmd)

rmse_results<-bind_rows(rmse_list_with_diff)
nlmd_reslts<-bind_rows(nlmd_list_with_diff)

all_results<-bind_rows(rmse_results, nlmd_reslts)

#missing data plot
ggplot(subset(all_results, missing_imputed=='missing'), aes(x=missingness, y=relative_diff, color=WV))+geom_boxplot()+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                     axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 18))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+ylab('Relative deviation from true WV estimates (%)')+
  xlab('Fraction of missing data')

#comparing imputed data

ggplot(subset(all_results, missing_imputed=='imputed'), aes(x=missingness, y=relative_diff, color=WV))+geom_boxplot()+
  facet_wrap(.~imputation_strategy, scales = "free")+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                     axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 18))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+ylab('Relative deviation from true WV estimates (%)')+
  xlab('Fraction of missing data')

#get means tables
WV_estimate_deviation_summaries_all<-all_results%>%
  group_by(WV, missing_imputed, missingness, imputation_strategy)%>%
  dplyr::summarise(mean_deviation=mean(relative_diff), se=se(relative_diff))%>%
  arrange(desc(missing_imputed), imputation_strategy, missingness, WV)

write.csv(WV_estimate_deviation_summaries_all, 'WV_estimation_summaries.csv', row.names = F)


################################
# true participant descriptives #
################################

#load df with true weights
descriptives<-readRDS("C:/Users/jaket/Dropbox/PhD/NoHoW Analyses/Weight variability/merged_test_2019-07-01/R_format/merged_nohow_dataset.RDS")%>%
  subset(., select=c(participant_study_id, elig_gender, elig_age, ecid1_weight_recorded, ecid1_height_recorded))%>%
  mutate(height_m=ecid1_height_recorded/100, BMI=ecid1_weight_recorded/(height_m)^2)%>%
  subset(select=c(-ecid1_height_recorded, -height_m))

ppt_info_included<-descriptives%>%subset(participant_study_id %in% df_true$ID)
colnames(ppt_info_included)<-c('ID', 'Gender', 'Age', 'Weight', 'BMI')
ppt_info_included[3:5]<-apply(ppt_info_included[3:5], 2, as.numeric)

ppt_info_included2<-df_true%>%group_by(ID)%>%
  filter(!is.na(weight))%>%dplyr::summarise(total_weights=n())%>%
  merge(., ppt_info_included, by='ID')

descrip_tab_all<-CreateTableOne(vars = c('Gender', 'Age', 'Weight', 'BMI', 'total_weights'),
                                data=ppt_info_included2)
print(descrip_tab_all)
 write.csv(print(descrip_tab_all), file="Protocol_paper_descriptives.csv")



