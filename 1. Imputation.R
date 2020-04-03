library(ggplot2)
library(dplyr)
library(lubridate)
library("birk")
library("imputeMulti")
library("imputeTS")
library(zoo)
library(xts)
library(reshape2)
library(plyr)
library(forecast)
library("tidyr")
library('rlist')
library(psych)
library('MLmetrics')
library('StatMeasures')
library('data.table')
library(anytime)
library(Metrics)
library('hydroGOF')
library(rowr)
library(simsem)
library(mice)
library(missForest)
library(gridExtra)

#load df
setwd("~/Dropbox/PhD/NoHoW Analyses/Weight variability")
#setwd("C:/Users/jaket/Dropbox/PhD/NoHoW Analyses/Weight variability")
df<-read.csv('most_complete5.csv')%>%arrange(ID, day_no) # read in from here

set.seed(123)
#function to insert random missingness in increments of 20, 40, 60 and 80%.
NAs_list<-list()
NA_generation_jitter<-function(x){
  for (i in 1:20){
    na20<-imposeMissing(data.mat=x, pmMCAR=0.2, ignoreCols=c('ID'))
    names(na20)[ncol(na20)] <- paste0("20_", i)
    na20<-na20[2]
    na40<-imposeMissing(data.mat=x, pmMCAR=0.4, ignoreCols=c('ID'))
    names(na40)[ncol(na40)] <- paste0("40_", i)
    na40<-na40[2]
    na60<-imposeMissing(data.mat=x, pmMCAR=0.6, ignoreCols=c('ID'))
    names(na60)[ncol(na60)] <- paste0("60_", i)
    na60<-na60[2]
    na80<-imposeMissing(data.mat=x, pmMCAR=0.8, ignoreCols=c('ID'))
    names(na80)[ncol(na80)] <- paste0("80_", i)
    na80<-na80[2]
    NAs_df<-cbind.data.frame(na20,na40,na60,na80)
    NAs_list[[i]]<-NAs_df
  }
  NAs_df<-cbind.data.frame(NAs_list)
  return(NAs_df)
} # This has created 4 (NA_20-80) * 20 (sims) * 50 (participants) simulated data sets


ID_weight<-subset(df, select=c(ID,weight)) #get minimal vars needed
set.seed(123)
NAs_list<-dlply(ID_weight, "ID", NA_generation_jitter)
NAs_df_MCAR<-plyr::ldply(NAs_list, cbind) #bind to get MCAR columns
NAs_list<-NULL

#write csv with MCAR simulated df
write.csv(NAs_df_MCAR, 'NA_df_50_MCAR_161219.csv', row.names = F)

#count the average NAs per each simulated column, we need this for plotting later
NA_count_total<-lapply(NAs_df_MCAR, function(x) sum(is.na(x)))%>%unlist()
NA_count_average<-data.frame(round(NA_count_total/50))
colnames(NA_count_average)<-'NAs'
NA_count_average$sim<-sub(".*_", "", rownames(NA_count_average))
NA_count_average$missingness<-as.factor(ifelse(grepl("80",rownames(NA_count_average)), '80',
                                               ifelse(grepl("60",rownames(NA_count_average)), '60',
                                                      ifelse(grepl("40",rownames(NA_count_average)), '40',
                                               ifelse(grepl("20",rownames(NA_count_average)), '20','0')))))
NA_count_average<-NA_count_average[-1,]
#write.csv(NA_count_average, 'NA_count_MCAR_291019.csv', row.names=F)

#check that the distribution of NAs looks okay (i.e. linear increases)
ggplot(NA_count_average, aes(x=missingness, y=NAs))+geom_point() # show dist of NAs

#sort df needed for next script on WV calculated
MCAR_df<-cbind.data.frame(df, NAs_df_MCAR)
MCAR_df[5]<-NULL


####################
# basic imputation #
####################
imputation_fn_big<-function(x){
  ts<-ts(x[,c(2:81)])
  tsclean_imp<-(apply (ts, 2, function(x) as.numeric(tsclean(x))))
  colnames(tsclean_imp)<-paste("TsC", colnames(tsclean_imp), sep="_")
  lin_int<-(apply (ts, 2, function(x) as.numeric(na.interpolation(x, option = "linear"))))
  colnames(lin_int)<-paste("Lint", colnames(lin_int), sep="_")
  stine_int<-(apply (ts, 2, function(x) as.numeric(na.interpolation(x, option = "stine"))))
  colnames(stine_int)<-paste("Stine", colnames(stine_int), sep="_")
  spline_int<-(apply (ts, 2, function(x) as.numeric(na.interpolation(x, option = "spline"))))
  colnames(spline_int)<-paste("Spline", colnames(spline_int), sep="_")
  kalman_impStr<-(apply (ts, 2, function(x) as.numeric(na.kalman(x, type="level", model = "StructTS"))))
  colnames(kalman_impStr)<-paste("StrKal", colnames(kalman_impStr), sep="_")
  kalman_impAuto<-(apply (ts, 2, function(x) as.numeric(na.kalman(x, model= 'auto.arima'))))
  colnames(kalman_impAuto)<-paste("AutoKal", colnames(kalman_impAuto), sep="_")
  EWMA<-(apply (ts, 2, function(x) as.numeric(na.ma(x, weighting  = "exponential", k=7))))
  colnames(EWMA)<-paste("EWMA", colnames(EWMA), sep="_")
  basic_imputations<-cbind.data.frame(x, tsclean_imp, lin_int, stine_int, spline_int, 
                                      kalman_impStr, kalman_impAuto, EWMA)
  return(basic_imputations)
}

#impuate all univariate NA columns
  MCAR_uni_imp<-dlply(MCAR_df, 'ID', imputation_fn_big)

#remove NA columns, leaving us only with imputed columns
  MCAR_uni_imp2<-lapply(MCAR_uni_imp, function(x) x=x[,c(1,82:641)])

#bind the cols for 1 big imputed df (removed NA columns)
  MCAR_uni_imp3<-plyr::ldply(MCAR_uni_imp2, rbind)
  MCAR_uni_imp3[1]<-NULL
  write.csv(MCAR_uni_imp3, 'basic_imp_big.csv', row.names=F) #write this so we dont have to re-do the imputation over

MCAR_uni_imp3<-read.csv('basic_imp_big.csv') 

## collect and merge real weights in
  MCAR_uni_imp4<-cbind(MCAR_uni_imp3, df)
  MCAR_uni_imp5<-MCAR_uni_imp4[,-c(1:2)]

#remove rows where weight was never there in true weight
  MCAR_uni_imp6<-filter(MCAR_uni_imp5, !is.na(weight))

# calculate the error between real and imputed dfs
real_rmse_list<-list()
rmse_fn<-function(x){
     for (i in 1:560){
         imputed<-as.vector(x[i])
         colnames(imputed)<-'not_real'
         real<-na.omit(x$weight)
         rmse<-Metrics::rmse(imputed$not_real, real)
         real_rmse_list[[i]]<-rmse
       }
     rmse_df<-bind_cols(real_rmse_list)
     return(rmse_df)
   }

#calculate RMSE between real and imputed dfs
mcar_uni_rmse_df<-dlply(MCAR_uni_imp6, 'ID', rmse_fn)%>%bind_rows()
colnames(mcar_uni_rmse_df)<-names1<-names(MCAR_uni_imp6)[1:560]

# grouping variables: amputation level (missingness) and imputation (strategy)
mcar_uni_rmse_df$missingness<-NA
mcar_uni_rmse_df$strategy<-NA

#df to long
mcar_uni_rmse_df_melt<-melt(mcar_uni_rmse_df)
colnames(mcar_uni_rmse_df_melt)[3:4]<-c('dataset', 'RMSE')

#assing missingness
mcar_uni_rmse_df_melt$missingness<-ifelse(grepl("_20_",mcar_uni_rmse_df_melt$dataset), 20, 
                                        ifelse(grepl("_40_",mcar_uni_rmse_df_melt$dataset), 40,
                                               ifelse(grepl("_60_",mcar_uni_rmse_df_melt$dataset), 60,80)))

mcar_uni_rmse_df_melt$missingness<-as.factor(mcar_uni_rmse_df_melt$missingness)

#assign strategy
mcar_uni_rmse_df_melt$strategy<-as.factor(ifelse(grepl("TsC_", mcar_uni_rmse_df_melt$dataset), 'TS Clean',
                               ifelse(grepl("Lint_", mcar_uni_rmse_df_melt$dataset), 'Linear Interpolation',
                                      ifelse(grepl("Stine_", mcar_uni_rmse_df_melt$dataset), 'Stine Interpolation',
                                             ifelse(grepl("Spline_", mcar_uni_rmse_df_melt$dataset), 'Spline Interpolation',
                                                    ifelse(grepl("StrKal_", mcar_uni_rmse_df_melt$dataset), 'StrucModel Kal',
                                                           ifelse(grepl("AutoKal_", mcar_uni_rmse_df_melt$dataset), 'Arima Kal',
                                                                         ifelse(grepl("EWMA_", mcar_uni_rmse_df_melt$dataset), 'Exp Weighted MA',
                                                                                'SeaDecom Kal'))))))))



##assign simulation number
mcar_uni_rmse_df_melt$sim<-sub(".*_", "", mcar_uni_rmse_df_melt$dataset)

##merge the real NA count in for plotting purposes
mcar_uni_rmse_df_melt2<-merge(mcar_uni_rmse_df_melt, NA_count_average, by=c('missingness', 'sim'))

#below we plot RMSE for MCAR univariate imps
gg_mcar_uni<-ggplot(data=mcar_uni_rmse_df_melt2, 
       aes(x=missingness, y=RMSE, color=missingness))+
  geom_boxplot()+facet_grid(.~strategy)+ylim(0,4)+
  theme(axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                 axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 14))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+ theme(legend.position = "none")

#take average errors per strategy and missingness
means_uni_mcar1<-mcar_uni_rmse_df_melt2%>%group_by(strategy, missingness)%>%
  dplyr::summarise(error=mean(RMSE), se=birk::se(RMSE))
means_uni_mcar1$insertion<-'MCAR'







####################################
####      REAL NA INSERTION     ###
###################################
# here we insert NA values in accordance with true patterns of missingness observed in real 
# participant's body weights. This improves external validity.

#get full data so we can find participants with increments of missingness
RPM_df<-read.csv("Daily_data_with_NAs_020719.csv")%>%filter(day_no<366)

#find n weights in 1 year, calculate it as a % of total missingness
percentage_na_all<-RPM_df%>%group_by(ID)%>%dplyr::summarise(weight_NAs=sum(is.na(weight)), total_days=n(),
                                                        perc.na=round((weight_NAs/total_days)*100),2)

#categorise these ppts into rough increments
percentage_na_all$group<-as.factor(ifelse(percentage_na_all$perc.na>78 & percentage_na_all$perc.na <82, 80,
                                ifelse(percentage_na_all$perc.na>58 & percentage_na_all$perc.na <62, 60,
                                       ifelse(percentage_na_all$perc.na>38 & percentage_na_all$perc.na <42, 40,
                                              ifelse(percentage_na_all$perc.na>18 & percentage_na_all$perc.na <22, 20, 0)))))

#select participants who only have missing data which aligns with necessary increments
percentage_na_interest<-filter(percentage_na_all, group !=0 & total_days>364)
percentage_na_interest$group<-as.factor(percentage_na_interest$group)

#take a random sample of ppts with eligible data missingness for each increment

random_sample_rpm_template_fn<-function(x){
  set.seed(123)
  samples<-x[sample(nrow(x), 20), ]
  samples$group<-x$group[1]
  return(samples)
}

rpm20<-dlply(percentage_na_interest, 'group', random_sample_rpm_template_fn)%>%bind_rows()

  write.csv(rpm20, 'rpm20.csv', row.names=F)
rpm20<-read.csv('rpm20.csv')

#we need to label each of these, quick fn to do it
add_no_missing_group_fn<-function(x){
  x$sim<-seq(1:nrow(x))
  return(x)
}
rpm20<-dlply(rpm20, 'group', add_no_missing_group_fn)%>%bind_rows()
IDs_for_NA<-levels(rpm20$ID) # collect their IDs

# now lets collect thereal data for these 20
real_data_for_RPM<-subset(RPM_df, ID %in% IDs_for_NA)%>%subset(., select=c(ID, day_no, weight))%>%merge(., rpm20, by='ID')%>%
  subset(., select=c(ID, day_no, group, sim, weight))%>%arrange(group, sim, ID)

#we want a col per ID, make binary for missingness
cast_fn<-function(x){
  out<-data.frame(x$weight)
  out$binary<-out[out > 0]<-1
  group<-as.character(x$group[1])
  sim<-as.character(x$sim[1])
  colnames(out)[1]<-paste(group, sim, sep="_")
  out[2]<-NULL
  return(out)
}
x<-filter(real_data_for_RPM, ID==real_data_for_RPM$ID[1])
rpm_list<-dlply(real_data_for_RPM, 'ID', cast_fn)
rpm_template<-cbind.data.frame(rpm_list)

#### we now have a template labelled. we need to impose each row on all FULL ppts.
# load in full ppts

#impose missingness patterns
imposing_list<-list()
imp_rmp_fn<-function(x){
  for (i in 1:80){
    imposing<-rpm_template
    imposing$ID<-x$ID[1]
    imposing[i]<-ifelse(imposing[i]==1, x$weight, NA)
    imposing_list[[i]]<-imposing[i]
  }
  imposed<-bind_cols(imposing_list)
  imposed$ID<-x$ID[1]
  imposed$day_no<-seq(1:nrow(imposed))
  return(imposed)
}

imposed_list<-dlply(df, 'ID', imp_rmp_fn)
imposed_df<-bind_rows(imposed_list)%>%dplyr::select(ID, day_no, everything())


#merge imposed with most complete so we have day no etc BY ID AND DAY NO
#all_rpm<-read.csv('~/Dropbox/PhD/NoHoW Analyses/Weight variability/most_complete5.csv')%>%
#  merge(., imposed_df, by=c('ID', 'day_no'))%>%arrange(ID, day_no)%>%subset(., select=c(-tseq))
all_rpm<-cbind(subset(df, select=weight), imposed_df)

# put an end row at the end of the ts of weights so we can impute.
add_end_row_fn<-function(x){
  ends<-data.frame(t(sapply(x, function(x) x[max(which(!is.na(x)))])))
  x[365,2:83]<-ends[,2:83]
  return(x)
}

all_rpm2<-dlply(all_rpm, 'ID', add_end_row_fn)%>%bind_rows()


#move ID, day_no etc to back so that 1:80 is the NA cols for imp
all_rpm3<-all_rpm[,c(4:83, 1:3)]

## how many NAs in each simulation??
NA_count_RPM<-sapply(all_rpm3[1:80], function(x) sum(is.na(x)))
NA_count_RPM<-data.frame(round(NA_count_RPM/50))
colnames(NA_count_RPM)[1]<-'NAs'
NA_count_RPM$sim<-sub(".*_", "", rownames(NA_count_RPM))
NA_count_RPM$missingness<-ifelse(grepl("20_",rownames(NA_count_RPM)), 20, 
                                 ifelse(grepl("40_",rownames(NA_count_RPM)), 40,
                                        ifelse(grepl("60_",rownames(NA_count_RPM)), 60,80)))

#ensure missingness % aligns with column names in RPM df
ggplot(NA_count_RPM, aes(x=missingness, y=NAs))+geom_point()


####################################
#     impute, univariate          #
###################################
#impute the RPM data using univariate algos
x<-filter(all_rpm3, ID==all_rpm3$ID[1])
uni_imp_fn2<-function(df){
  ID<-df$ID[1]
  df$ID<-NULL
  sim<-colnames(df)[1]
  colnames(df)[1]<-'imputed_weight'
  tsclean_imp<-data.frame(try(tsclean(df$imputed_weight)))
  lin_int<-data.frame((try(na.interpolation(df$imputed_weight, option = "linear"))))
  stine_int<-data.frame((try(na.interpolation(df$imputed_weight, option = "stine"))))
  spline_int<-data.frame((try(na.interpolation(df$imputed_weight, option = "spline"))))
  kalman_impStr<-data.frame((try(na.kalman(df$imputed_weight, type="level", model = "StructTS"))))
  kalman_impAuto<-data.frame((try(na.kalman(df$imputed_weight, model= 'auto.arima'))))
  EWMA<-data.frame((try(na.ma(df$imputed_weight, weighting  = "exponential", k=7))))
  basic_imputations<-bind_cols(tsclean_imp, lin_int, stine_int, spline_int, 
                               kalman_impStr, kalman_impAuto, EWMA)
  colnames(basic_imputations)<-c(paste('TsC', sim, sep="_"), paste("Lint", sim, sep="_"), paste("Stine", sim, sep="_"),
                                 paste("Spline", sim, sep="_"), paste("StrKal", sim, sep="_"), 
                                 paste("AutoKal", sim, sep="_"),paste("EWMA", sim, sep="_"))
  return(basic_imputations)
}


#for all IDs, impute all columns
out=list()
i=1
multi_loop_fn<-function(x){
  for (i in 1:80){
    df<-x[, c(i,81, 82)]
    results<-uni_imp_fn2(df)
    out[[i]]<-results
  }
  out_df<-bind_cols(out)
  return(out_df)
}


  RPM_univariate_imp_df<-dlply(all_rpm3, 'ID', multi_loop_fn)%>%bind_rows()

#write this so we dont have to repeat due to time
  write.csv(RPM_univariate_imp_df, 'RPM_univariate_imp_df.csv', row.names=F) 
  
RPM_univariate_imp_df<-read.csv('RPM_univariate_imp_df.csv') #then just read it in

RPM_univariate_imp_df2<-cbind(RPM_univariate_imp_df, subset(df, select=c(-tseq))) #add ID, day_no, weight back in
RPM_univariate_imp_df3<-filter(RPM_univariate_imp_df2, !is.na(weight))

#calculate RMSE
real_rmse_list<-list()
rmse_fn2<-function(x){
  for (i in 1:560){
    rmse<-as.numeric(hydroGOF::rmse(x[i], x$weight, na.rm=T))
    real_rmse_list[[i]]<-rmse
  }
  rmse_df<-bind_cols(real_rmse_list)
  return(rmse_df)
}

#sort the structure of the df to something less retarded
RMSE_RPM_uni_df<-dlply(RPM_univariate_imp_df3, 'ID', rmse_fn2)%>%bind_rows()
colnames(RMSE_RPM_uni_df)<-names(RPM_univariate_imp_df2)[1:560]

# grouping variables: amputation level (missingness) and imputation (strategy)
RMSE_RPM_uni_df$missingness<-NA
RMSE_RPM_uni_df$strategy<-NA

#df to long
RMSE_RPM_uni_df_melt<-melt(RMSE_RPM_uni_df)
colnames(RMSE_RPM_uni_df_melt)[3:4]<-c('dataset', 'RMSE')

#assing missingness
RMSE_RPM_uni_df_melt$missingness<-ifelse(grepl("_20_",RMSE_RPM_uni_df_melt$dataset), 20, 
                                         ifelse(grepl("_40_",RMSE_RPM_uni_df_melt$dataset), 40,
                                                ifelse(grepl("_60_",RMSE_RPM_uni_df_melt$dataset), 60,80)))

RMSE_RPM_uni_df_melt$missingness<-as.factor(RMSE_RPM_uni_df_melt$missingness)
#assign strategy
RMSE_RPM_uni_df_melt$strategy<-ifelse(grepl("TsC_", RMSE_RPM_uni_df_melt$dataset), 'TS Clean',
                                      ifelse(grepl("Lint_", RMSE_RPM_uni_df_melt$dataset), 'Linear Interpolation',
                                             ifelse(grepl("Stine_", RMSE_RPM_uni_df_melt$dataset), 'Stine Interpolation',
                                                    ifelse(grepl("Spline_", RMSE_RPM_uni_df_melt$dataset), 'Spline Interpolation',
                                                           ifelse(grepl("StrKal_", RMSE_RPM_uni_df_melt$dataset), 'StrucModel Kal',
                                                                  ifelse(grepl("AutoKal_", RMSE_RPM_uni_df_melt$dataset), 'Arima Kal',
                                                                         ifelse(grepl("EWMA_", RMSE_RPM_uni_df_melt$dataset), 'Exp Weighted MA',
                                                                                NA)))))))
RMSE_RPM_uni_df_melt$strategy<-as.factor(RMSE_RPM_uni_df_melt$strategy)
##assign simulation number
RMSE_RPM_uni_df_melt$sim<-sub(".*_", "", RMSE_RPM_uni_df_melt$dataset)

#######################
# merge into NA count #
#######################

RMSE_RPM_uni_df_melt2<-merge(RMSE_RPM_uni_df_melt, NA_count_RPM, by=c('missingness', 'sim'))

t<-RMSE_RPM_uni_df_melt2%>%group_by(missingness)%>%dplyr::summarise(mean=mean(RMSE))

gg_rpm_uni<-ggplot(data=RMSE_RPM_uni_df_melt2, aes(x=missingness, y=RMSE, color=missingness))+
  geom_boxplot()+facet_grid(.~strategy)+ylim(0,4)+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                     axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 14))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+xlab('Number of missing values')+ theme(legend.position = "none")



means_rpm_uni1<-RMSE_RPM_uni_df_melt2%>%group_by(strategy, missingness)%>%
  dplyr::summarise(error=mean(RMSE), se=birk::se(RMSE))
means_rpm_uni1$insertion<-'RPM'


uni_imp_RMSE<-bind_rows(means_rpm_uni1, means_uni_mcar1)
write.csv(uni_imp_RMSE, 'rpm_uni_all.csv', row.names = F)











###############################
# multivariate imputation #####
##############################
#
multi_mcar_df<-MCAR_df%>%mutate(tseq=anydate(tseq), dow=weekdays(tseq))%>%subset(select=c(-tseq, -weight))
multi_mcar_df2<-multi_mcar_df[,c(3:82, 1:2, 83)]

x<-filter(multi_mcar_df2, ID==multi_mcar_df2$ID[1])
#run 3 multivar imps for each na col
multi_imp_fn<-function(df){
  ID<-df$ID
  df$ID<-NULL
  df$dow<-as.factor(df$dow)
  sim<-colnames(df)[1]
  colnames(df)[1]<-'weight'
  imp_KNN<-data.frame(DMwR::knnImputation(df, k=2)[1])
  imp_mice_pmm<-data.frame(mice(df, m=10)%>%mice::complete(., action='all')%>%
                             bind_rows(.)%>%group_by(day_no)%>%dplyr::summarise(weight=mean(weight))%>%
                             subset(., select=c(-day_no)))
  imp_RF<-data.frame(missForest(df, ntree=100)$ximp[1])
  multi_var_imp<-cbind(imp_KNN, imp_mice_pmm, imp_RF)
  colnames(multi_var_imp)<-c(paste('KNN', sim, sep='_'), paste('PMM', sim, sep='_'), paste('RF', sim, sep='_'))
  return(multi_var_imp)
}

multi_loop_fn<-function(x){
  for (i in 1:80){
    df<-x[, c(i, 81:83)]
    results<-multi_imp_fn(df)
    out[[i]]<-results
  }
  out_df<-bind_cols(out)
  return(out_df)
}

  multi_imp_df<-dlply(multi_mcar_df2, 'ID', multi_loop_fn)%>%bind_rows()
  write.csv(multi_imp_df, 'multi_imp_MAPE_df.csv', row.names = F)

multi_imp_df<-read.csv('multi_imp_MAPE_df.csv')
  
#merge with observed data
multi_imp_df2<-bind_cols(multi_imp_df, df)
  
#calculate RMSE
real_rmse_list<-list()
rmse_fn3<-function(x){
  for (i in 1:240){
    rmse<-as.numeric(hydroGOF::rmse(x[i], x$weight, na.rm=T))
    real_rmse_list[[i]]<-rmse
  }
  rmse_df<-bind_cols(real_rmse_list)
  return(rmse_df)
}

rmse_multi_mcar_df<-dlply(multi_imp_df2, 'ID', rmse_fn3)%>%bind_rows()
colnames(rmse_multi_mcar_df)<-colnames(multi_imp_df2)[1:240]


# grouping variables: amputation level (missingness) and imputation (strategy)
rmse_multi_mcar_df$missingness<-NA
rmse_multi_mcar_df$strategy<-NA

#df to long
rmse_multi_mcar_melt<-melt(rmse_multi_mcar_df)
colnames(rmse_multi_mcar_melt)[3:4]<-c('dataset', 'RMSE')

#assing missingness
rmse_multi_mcar_melt$missingness<-as.factor(ifelse(grepl("_20_",rmse_multi_mcar_melt$dataset), 20, 
                                 ifelse(grepl("_40_",rmse_multi_mcar_melt$dataset), 40,
                                        ifelse(grepl("_60_",rmse_multi_mcar_melt$dataset), 60,80))))

#assign strategy
rmse_multi_mcar_melt$strategy<-as.factor(ifelse(grepl("KNN",rmse_multi_mcar_melt$dataset), 'KNN', 
                              ifelse(grepl("PMM",rmse_multi_mcar_melt$dataset), 'PMM',
                                     ifelse(grepl("RF",rmse_multi_mcar_melt$dataset), 'RF','sim'))))


##assign simulation number
rmse_multi_mcar_melt$sim<-sub(".*_", "", rmse_multi_mcar_melt$dataset)

##merge the real NA count in for plotting purposes
rmse_multi_mcar_melt2<-merge(rmse_multi_mcar_melt, NA_count_average, by=c('missingness', 'sim'))


#below we plot RMSE for MCAR univariate imps
gg_mcar_multi<-ggplot(data=rmse_multi_mcar_melt2, 
       aes(x=missingness, y=RMSE, color=missingness))+
  geom_boxplot()+facet_grid(.~strategy)+ylim(0,4)+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                     axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 18))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+ theme(legend.position = "none")

#take average errors per strategy and missingness
means_multi_mcar1<-rmse_multi_mcar_melt2%>%group_by(strategy, missingness)%>%
  dplyr::summarise(error=mean(RMSE), se=birk::se(RMSE))
means_multi_mcar1$insertion<-'MCAR'




########################
# multi RPM imputation #
########################

#sort RPM df for multi imputation
multi_rpm_na_df<-all_rpm3[,-c(81:83)]
multi_rpm_na_df2<-cbind(multi_rpm_na_df, df)%>%mutate(tseq=anydate(tseq), dow=weekdays(tseq))%>%subset(select=c(-tseq, -weight))

#multi impute
  multi_imp_rpm_df<-dlply(multi_rpm_na_df2, 'ID', multi_loop_fn)%>%bind_rows()
  
write.csv(multi_imp_rpm_df, 'multi_imp_rpm_df.csv', row.names = F)
multi_imp_rpm_df<-read.csv('multi_imp_rpm_df.csv')

#merge with observed data
multi_imp_rpm_df2<-bind_cols(multi_imp_rpm_df, df)

#calculate RMSE
real_rmse_list<-list()
rmse_multi_rpm_df<-dlply(multi_imp_rpm_df2, 'ID', rmse_fn3)%>%bind_rows()
colnames(rmse_multi_rpm_df)<-colnames(multi_imp_rpm_df2)[1:240]

# grouping variables: amputation level (missingness) and imputation (strategy)
rmse_multi_rpm_df$missingness<-NA
rmse_multi_rpm_df$strategy<-NA

#df to long
rmse_multi_rpm_melt<-melt(rmse_multi_rpm_df)
colnames(rmse_multi_rpm_melt)[3:4]<-c('dataset', 'RMSE')

#assing missingness
rmse_multi_rpm_melt$missingness<-as.factor(ifelse(grepl("_20_",rmse_multi_rpm_melt$dataset), 20, 
                                                   ifelse(grepl("_40_",rmse_multi_rpm_melt$dataset), 40,
                                                          ifelse(grepl("_60_",rmse_multi_rpm_melt$dataset), 60,80))))

#assign strategy
rmse_multi_rpm_melt$strategy<-as.factor(ifelse(grepl("KNN",rmse_multi_rpm_melt$dataset), 'KNN', 
                                                ifelse(grepl("PMM",rmse_multi_rpm_melt$dataset), 'PMM',
                                                       ifelse(grepl("RF",rmse_multi_rpm_melt$dataset), 'RF','sim'))))


##assign simulation number
rmse_multi_rpm_melt$sim<-sub(".*_", "", rmse_multi_rpm_melt$dataset)

##merge the real NA count in for plotting purposes
rmse_multi_rpm_melt2<-merge(rmse_multi_rpm_melt, NA_count_average, by=c('missingness', 'sim'))


#below we plot RMSE for MCAR univariate imps
gg_rpm_multi<-ggplot(data=rmse_multi_rpm_melt2, 
       aes(x=missingness, y=RMSE, color=missingness))+
  geom_boxplot()+facet_grid(.~strategy)+ylim(0,4)+
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+ theme(axis.text.x=element_text(size=16,color='black', angle=90),
                                                     axis.text.y=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 18))+ 
  theme(legend.title = element_text(color = "black", size = 18),
        legend.text = element_text(color = "black"))+ 
  theme(legend.text = element_text(colour="black", size=18))+ theme(legend.position = "none")

#take average errors per strategy and missingness
means_multi_rpm1<-rmse_multi_rpm_melt2%>%group_by(strategy, missingness)%>%
  dplyr::summarise(error=mean(RMSE), se=birk::se(RMSE))
means_multi_rpm1$insertion<-'RPM'


#list of main imputation plots
gg_list<-list(gg_mcar_uni, gg_rpm_uni, gg_mcar_multi, gg_rpm_multi)

do.call(grid.arrange,gg_list) ## facet these 4 grips from list (::gridExtra)



#aggregate RMSE values for comparative imputation table
summary_rmse<-bind_rows(means_multi_rpm1, means_uni_mcar1, means_rpm_uni1, means_multi_mcar1)
summary_rmse[3:4]<-round(summary_rmse[3:4], 3)

write.csv(summary_rmse, 'imputation_rmse_summary_table.csv', row.names=F)






##########################
# example plot for paper #
##########################
levels(MCAR_uni_imp6$ID)
ppt<-filter(multi_imp_df, ID=='CPH_EMN_7406')
ppt$strategy<-NA
ppt$missingness<-NA
ppt$ID<-NULL

ppt_melt<-melt(ppt, id.vars=c('strategy', 'missingness', 'day_no', 'tseq'))

ppt_melt$strategy<-ifelse(grepl("TsC_", ppt_melt$variable), 'TS Clean',
                              ifelse(grepl("Lint_", ppt_melt$variable), 'Linear Interpolation',
                                     ifelse(grepl("Stine_", ppt_melt$variable), 'Stine Interpolation',
                                            ifelse(grepl("Spline_", ppt_melt$variable), 'Spline Interpolation',
                                                   ifelse(grepl("StrKal_", ppt_melt$variable), 'StrucModel Kal',
                                                          ifelse(grepl("AutoKal_", ppt_melt$variable), 'Arima Kal',
                                                                 ifelse(grepl("EWMA_", ppt_melt$variable), 'Exp Weighted MA',
                                                                        ifelse(grepl("weight", ppt_melt$variable), 'Weight',
                                                                               'SeaDecom Kal'))))))))
ppt_melt$missingness<-ifelse(grepl("_20",ppt_melt$variable), 20, 
                                 ifelse(grepl("_40",ppt_melt$variable), 40,
                                        ifelse(grepl("_60",ppt_melt$variable), 60,
                                               ifelse(grepl("_80",ppt_melt$variable), 80, 0))))


### Structural model - kalman smoother plots ##

ppt_melt%>%filter(variable=='weight' | variable=='StrKal_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("SMKS 40%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('Structural Model Kalman Smoother - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='StrKal_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("SMKS 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Structural Model Kalman Smoother - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## EWMA plots ##

ppt_melt%>%filter(variable=='weight' | variable=='EWMA_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("EWMA 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Exponentially Weighted Moving Average - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='EWMA_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("EWMA 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Exponentially Weighted Moving Average - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## linear interpolation ##

ppt_melt%>%filter(variable=='weight' | variable=='Lint_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Lin_inter 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Linear interpolation - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='Lint_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Lin_inter 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Linear interpolation - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## stineman interpolation ##

ppt_melt%>%filter(variable=='weight' | variable=='Stine_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Stine_inter 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Stineman interpolation - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='Stine_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Stine_inter 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Stineman interpolation - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## TSClean ##

ppt_melt%>%filter(variable=='weight' | variable=='TsC_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Ts Clean 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Ts Clean - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))


ppt_melt%>%filter(variable=='weight' | variable=='TsC_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Ts Clean 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Ts Clean - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## ARIMA Kalman ##

ppt_melt%>%filter(variable=='weight' | variable=='AutoKal_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("ARIMA_Kal 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('ARIMA state-space representation and Kalman smoothing - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='AutoKal_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("ARIMA_Kal 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('ARIMA state-space representation and Kalman smoothing - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## spline interpolation ##

ppt_melt%>%filter(variable=='weight' | variable=='Spline_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Spline_int 40%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Spline interpolation - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_melt%>%filter(variable=='weight' | variable=='Spline_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("Spline_int 80%", "Real" ), values=c('tomato','steelblue1'))+
  ggtitle('Spline interpolation - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

################################
## multivariate insertion

ppt_multi<-filter(multi_imp_df2, ID=='CPH_EMN_7406')
ppt_multi$strategy<-NA
ppt_multi$missingness<-NA
ppt_multi$ID<-NULL

ppt_multi_melt<-melt(ppt_multi, id.vars=c('strategy', 'missingness', 'day_no', 'tseq'))

ppt_multi_melt$strategy<-as.factor(ifelse(grepl("KNN",ppt_multi_melt$variable), 'KNN', 
                                            ifelse(grepl("PMM",ppt_multi_melt$variable), 'PMM',
                                                   ifelse(grepl("RF",ppt_multi_melt$variable), 'RF','sim'))))

ppt_multi_melt$missingness<-ifelse(grepl("_20",ppt_multi_melt$variable), 20, 
                             ifelse(grepl("_40",ppt_multi_melt$variable), 40,
                                    ifelse(grepl("_60",ppt_multi_melt$variable), 60,
                                           ifelse(grepl("_80",ppt_multi_melt$variable), 80, 0))))


#PMM##

ppt_multi_melt%>%filter(variable=='weight' | variable=='PMM_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("PMM 40%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('Predictive Means Matching - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_multi_melt%>%filter(variable=='weight' | variable=='PMM_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("PMM 80%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('Predictive Means Matching - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## KNN ##

ppt_multi_melt%>%filter(variable=='weight' | variable=='KNN_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("KNN 40%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('K-Nearest Neighbour - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_multi_melt%>%filter(variable=='weight' | variable=='KNN_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("KNN 80%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('K-Nearest Neighbour - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

## RF ##

ppt_multi_melt%>%filter(variable=='weight' | variable=='RF_40_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("RF 40%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('Random forest - 40% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

ppt_multi_melt%>%filter(variable=='weight' | variable=='RF_80_1')%>% 
  ggplot(., aes(y=value, x=as.Date(tseq)))+geom_point()+geom_line(aes(colour=strategy), size=1)+ 
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16))+ 
  theme(axis.text=element_text(size=16))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(legend.title = element_blank(), legend.text = element_text(size=16))+xlab('Date')+ylab('Weight(kg)')+
  scale_color_manual(values=c('tomato','steelblue1'))+
  labs(color = "Legend Title\n") +
  scale_color_manual(labels = c("RF 80%", "Real"), values=c('tomato','steelblue1'))+
  ggtitle('Random forest - 80% imputed')+
  theme(plot.title = element_text(hjust = 0.5))

########################################
# plotting for simulation illustration #
########################################
#get ppt ppt
ppt<-filter(all_rpm, ID=='CPH_EMN_7406')
colnames(ppt)[8]<-'fourty_ex'
colnames(ppt)[14]<-'eighty_ex'

#observed weight
ggplot(ppt, aes(x=day_no, y=weight))+geom_point()

#40% missing
ggplot(ppt, aes(x=day_no, y=fourty_ex))+geom_point()

#80% missing
ggplot(ppt, aes(x=day_no, y=eighty_ex))+geom_point()+

