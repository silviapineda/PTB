rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB Diags
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Analysis of EMR data for PTB
###         
###
### Author: Silvia Pineda
### Date: April, 2018
############################################################################################
library(lattice)
library(lme4)
library("RColorBrewer")
library(ggplot2)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)


EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")


#########################
#### Diags data ########
########################
##Extract Individual ID for the Patient_index
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_diags$Patient_index),"_")
EMR_long_diags$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_diags$Birth_id<- unlist(lapply(splitpop, "[", 2))
EMR_long_diags$Unique_id<-paste(EMR_long_diags$Patient_index,EMR_long_diags$WeekOfPregnancy,sep="_")

##Select all the variables 
EMR_long_diags_result<-EMR_long_diags[,7:(ncol(EMR_long_diags)-3)] #3157 diags
rownames(EMR_long_diags_result)<-EMR_long_diags$Unique_id
table(EMR_long_diags$Term) #1928 PTB and 32291 Term

matrix<-data.matrix(table(EMR_long_diags$Patient_index,EMR_long_diags$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #414 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #8648 births with Term

###Select the lab test that are complete
num_null_TERM<-NULL
num_null_PTB<-NULL
for (i in 1:ncol(EMR_long_diags_result)){
  num_null_TERM[i]<-table(EMR_long_diags_result[which(EMR_long_diags$Term=="Term"),i])[2]
num_null_PTB[i]<-table(EMR_long_diags_result[which(EMR_long_diags$Term=="PTB"),i])[2]
}

EMR_long_diags_result_filter_PTB<-EMR_long_diags_result[which(EMR_long_diags$Term=="PTB"),which(num_null_PTB>2)] #503 that has more than 2 samples with result
EMR_long_diags_result_filter_Term<-EMR_long_diags_result[which(EMR_long_diags$Term=="Term"),which(num_null_TERM>2)] #1576 that has more than 2 samples with result

id.match<-match(colnames(EMR_long_diags_result_filter_PTB),colnames(EMR_long_diags_result_filter_Term))
PTB_data<-EMR_long_diags_result_filter_PTB[,which(is.na(id.match)==F)]
Term_data<-EMR_long_diags_result_filter_Term[,na.omit(id.match)]
EMR_long_diags_result_filter<-rbind(PTB_data,Term_data)
id.merge<-match(EMR_long_diags$Unique_id,rownames(EMR_long_diags_result_filter))
EMR_long_diags_result_filter<-EMR_long_diags_result_filter[id.merge,]
##473 diags final set


####Running the longitudinal model
EMR_long_diags$Patient_index<-factor(EMR_long_diags$Patient_index)
EMR_long_diags$Term<-factor(EMR_long_diags$Term,levels = c("Term","PTB"))
EMR_long_diags$Individual_id<-factor(EMR_long_diags$Individual_id)
EMR_long_diags$Outcome<-ifelse(EMR_long_diags$Term=="PTB",1,0)

results_diags<-matrix(NA,ncol(EMR_long_diags_result_filter),2)
for (i in 1:ncol(EMR_long_diags_result_filter)){
  print(i)
  fm_full <-  try(glmer(EMR_long_diags$Outcome ~ EMR_long_diags_result_filter[,i] +
                          (1|EMR_long_diags$Patient_index),
                        family=binomial))
  
  if(class(fm_full)!="try-error"){
    results_diags[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
    results_diags[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
  }
}
colnames(results_diags)<-c("coef_diags","p_diags")
rownames(results_diags)<-colnames(EMR_long_diags_result_filter)
write.csv(results_diags,"results_diags.csv")

results_diags<-read.csv("Results/DIAGS/results_diags.csv")
padjusted<-p.adjust(results_diags[,3],"fdr")
results_diags_sign<-results_diags[which(padjusted<0.05),] #196 significance(all pass MT)
id.sign<-match(results_diags_sign$X,colnames(EMR_long_diags_result_filter))
EMR_long_diags_result_filter_sign<-EMR_long_diags_result_filter[,id.sign]

perc_term<-NULL
perc_PTB<-NULL
num_term<-NULL
num_PTB<-NULL
for (i in 1:ncol(EMR_long_diags_result_filter_sign)){
  num_term[i]<-table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,1]
  num_PTB[i]<-table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,2]
  perc_term[i]<-(table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,]/ table(EMR_long_diags$Term))[1]
  perc_PTB[i]<-(table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,]/ table(EMR_long_diags$Term))[2]
  test_Term<-EMR_long_diags_result_filter_sign[which(EMR_long_diags$Term=="Term"),i]
  test_PTB<-EMR_long_diags_result_filter_sign[which(EMR_long_diags$Term=="PTB"),i]
}

write.csv(cbind(results_diags_sign,num_term,perc_term,num_PTB,perc_PTB),file="Results/DIAGS/results_diags_sign.csv")

####To plot the data
for(i in 1:ncol(EMR_long_diags_result_filter_sign)){
  print(i)
  EMR_long_diags$diags<-EMR_long_diags_result_filter_sign[,i]
  EMR_long_diags$Term<-factor(EMR_long_diags$Term,levels = c("PTB","Term"))
  
  
  tiff(paste0("Results/DIAGS/",colnames(EMR_long_diags_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
  p2<-ggplot(EMR_long_diags, aes(x=as.character(WeekOfPregnancy))) +
    geom_bar(data=EMR_long_diags[EMR_long_diags$Term=="Term",],
             aes(y=(diags)/length(diags),fill=Term), stat="identity") +
    geom_bar(data=EMR_long_diags[EMR_long_diags$Term=="PTB",],
             aes(y=-(diags)/length(diags),fill=Term), stat="identity") +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip(ylim=c(-0.1,0.1)) +
    scale_y_continuous(breaks=seq(-0.1,0.1,0.05), labels=c(0.1,0.05,0,0.05,0.1)) +
    labs(y="Percentage of samples with diagnosis", x="Week of pregnancy") +
    ggtitle("                         PTB                                      Term")
  print(p2)
  dev.off()
}


#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
EMR_long_diags_multivariate<-cbind(EMR_long_diags$Outcome,EMR_long_diags$Patient_index,EMR_long_diags_result_filter_sign)
colnames(EMR_long_diags_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_diags_multivariate,file="EMR_long_diags_multi.Rdata")

##To extract the names of the variables
paste(colnames(EMR_long_diags_multivariate)[1:50],collapse ="+")
paste(colnames(EMR_long_diags_multivariate)[51:100],collapse ="+")
paste(colnames(EMR_long_diags_multivariate)[101:198],collapse ="+")

fm_full <-  glmer(Term ~ A15.9+A56.19+A60.04+A92.5+B06.9+B96.89+D18.01+D23.9+D63.1+D69.6+D83.9+E01.0+E05.90+E06.9+E10.65+E11.319+E11.39+E23.6+E66.3+E78.5+E87.5+E88.09+E89.0+F11.10+F12.20+F31.9+F43.10+F43.20+F53+F90.9+G40.209+G40.909+G43.009+G43.109+G44.219+G47.00+G51.0+G56.03+H20.9+H93.11+H93.19+I15.8+I37.0+I42.0+I42.4+I45.81+I51.9+J02.9+
                    J30.89+J32.9+J45.30+J98.11+K21.9+K31.89+K43.9+K51.90+K63.89+K90.0+L71.0+L71.9+L72.3+L73.8+L73.9+L81.9+M06.9+M25.531+M25.539+M25.569+M54.16+M54.2+M54.32+M62.89+M79.1+N05.8+N10+N18.2+N18.3+N20.0+N28.9+N30.10+N86+N97.0+O00.90+O09.02+O09.40+O09.73+O09.891+O10.419+O10.912+O11.9+O13.1+O14.93+O16.1+O22.33+O23.02+O23.03+O26.02+O26.10+
                    O26.12+O30.009+O30.049+O30.109+O31.10X0+O33.7XX0+O34.02+O34.29+O34.43+O35.2XX0+O36.0120+O36.0199+O36.5921+O36.5991+O41.03X1+O41.8X20+O41.8X30+O41.90X0+O43.129+O44.12+O69.0XX0+O90.3+O92.29+O92.3+O92.79+O98.312+O98.313+O9A.319+P01.3+P02.29+P03.819+Q05.9+Q20.1+Q21.1+Q28.3+Q33.0+Q51.818+Q62.5+Q89.9+R06.00+R06.09+R09.89+R10.31+R11.10+R18.8+R19.09+R19.8+R20.9+R22.0+R22.1+R31.9+R32+R58+R63.4+R63.5+R63.6+R71.8+R76.0+R76.8+R82.79+R91.8+S80.00XA+Y92.099+Z03.71+Z03.89+Z12.4+Z20.1+Z30.2+Z31.430+Z36.82+Z38.2+Z39.2+Z48.22+Z51.81+Z53.1+Z59.1+Z59.9+Z65.9+Z67.41+Z68.26+Z68.27+Z68.33+Z68.36+Z68.42+Z80.8+Z80.9+Z82.0+Z84.81+Z85.850+Z87.19+Z87.39+Z87.74+Z88.0+Z90.49+Z91.040+Z91.19+Z91.49+Z94.1+
                    (1|Patient_index),family=binomial,data = EMR_long_diags_multivariate)

#########################
##  Run in the server  ##

load("EMR_long_diags_multi.Rdata")
library(glmmLasso)

##Number of unique patient_index to obtain the train and test set
set.seed(33)
EMR_long_diags_patient<-EMR_long_diags_multivariate[unique(EMR_long_diags_multivariate$Patient_index),c("Term","Patient_index")]
table(EMR_long_diags_patient[,"Term"])
##20% of Term and PTB
predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:442),88)
  EMR_test_PTB<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:8620),1724)
  EMR_test_Term<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==0)[id_test_data_Term],]
  
  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_diags_multivariate$Patient_index)))
  EMR_long_diags_test<-EMR_long_diags_multivariate[id_test_data,] 
  EMR_long_diags_train<-EMR_long_diags_multivariate[-(id_test_data),] 
  
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  ## first fit good starting model
  library(MASS);library(nlme)
  #PQL<-glmmPQL(Term~1,random = ~1|Patient_index,family=family,data=EMR_long_diags_train)
  
  #H93.19
  EMR_long_diags_train<-EMR_long_diags_train[,-43]
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ A15.9+A56.19+A60.04+A92.5+B06.9+B96.89+D18.01+D23.9+D63.1+D69.6+D83.9+E01.0+E05.90+E06.9+E10.65+E11.319+E11.39+E23.6+E66.3+E78.5+E87.5+E88.09+E89.0+F11.10+F12.20+F31.9+F43.10+F43.20+F53+F90.9+G40.209+G40.909+G43.009+G43.109+G44.219+G47.00+G51.0+G56.03+H20.9+H93.11+I15.8+I37.0+I42.0+I42.4+I45.81+I51.9+J02.9+
                            J30.89+J32.9+J45.30+J98.11+K21.9+K31.89+K43.9+K51.90+K63.89+K90.0+L71.0+L71.9+L72.3+L73.8+L73.9+L81.9+M06.9+M25.531+M25.539+M25.569+M54.16+M54.2+M54.32+M62.89+M79.1+N05.8+N10+N18.2+N18.3+N20.0+N28.9+N30.10+N86+N97.0+O00.90+O09.02+O09.40+O09.73+O09.891+O10.419+O10.912+O11.9+O13.1+O14.93+O16.1+O22.33+O23.02+O23.03+O26.02+O26.10+
                            O26.12+O30.009+O30.049+O30.109+O31.10X0+O33.7XX0+O34.02+O34.29+O34.43+O35.2XX0+O36.0120+O36.0199+O36.5921+O36.5991+O41.03X1+O41.8X20+O41.8X30+O41.90X0+O43.129+O44.12+O69.0XX0+O90.3+O92.29+O92.3+O92.79+O98.312+O98.313+O9A.319+P01.3+P02.29+P03.819+Q05.9+Q20.1+Q21.1+Q28.3+Q33.0+Q51.818+Q62.5+Q89.9+R06.00+R06.09+R09.89+R10.31+R11.10+R18.8+R19.09+R19.8+R20.9+R22.0+R22.1+R31.9+R32+R58+R63.4+R63.5+R63.6+R71.8+R76.0+R76.8+R82.79+R91.8+S80.00XA+Y92.099+Z03.71+Z03.89+Z12.4+Z20.1+Z30.2+Z31.430+Z36.82+Z38.2+Z39.2+Z48.22+Z51.81+Z53.1+Z59.1+Z59.9+Z65.9+Z67.41+Z68.26+Z68.27+Z68.33+Z68.36+Z68.42+Z80.8+Z80.9+Z82.0+Z84.81+Z85.850+Z87.19+Z87.39+Z87.74+Z88.0+Z90.49+Z91.040+Z91.19+Z91.49+Z94.1,
                            rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  print(opt)
  glm_final <- try(glmmLasso(Term ~ A15.9+A56.19+A60.04+A92.5+B06.9+B96.89+D18.01+D23.9+D63.1+D69.6+D83.9+E01.0+E05.90+E06.9+E10.65+E11.319+E11.39+E23.6+E66.3+E78.5+E87.5+E88.09+E89.0+F11.10+F12.20+F31.9+F43.10+F43.20+F53+F90.9+G40.209+G40.909+G43.009+G43.109+G44.219+G47.00+G51.0+G56.03+H20.9+H93.11+I15.8+I37.0+I42.0+I42.4+I45.81+I51.9+J02.9+
                           J30.89+J32.9+J45.30+J98.11+K21.9+K31.89+K43.9+K51.90+K63.89+K90.0+L71.0+L71.9+L72.3+L73.8+L73.9+L81.9+M06.9+M25.531+M25.539+M25.569+M54.16+M54.2+M54.32+M62.89+M79.1+N05.8+N10+N18.2+N18.3+N20.0+N28.9+N30.10+N86+N97.0+O00.90+O09.02+O09.40+O09.73+O09.891+O10.419+O10.912+O11.9+O13.1+O14.93+O16.1+O22.33+O23.02+O23.03+O26.02+O26.10+
                           O26.12+O30.009+O30.049+O30.109+O31.10X0+O33.7XX0+O34.02+O34.29+O34.43+O35.2XX0+O36.0120+O36.0199+O36.5921+O36.5991+O41.03X1+O41.8X20+O41.8X30+O41.90X0+O43.129+O44.12+O69.0XX0+O90.3+O92.29+O92.3+O92.79+O98.312+O98.313+O9A.319+P01.3+P02.29+P03.819+Q05.9+Q20.1+Q21.1+Q28.3+Q33.0+Q51.818+Q62.5+Q89.9+R06.00+R06.09+R09.89+R10.31+R11.10+R18.8+R19.09+R19.8+R20.9+R22.0+R22.1+R31.9+R32+R58+R63.4+R63.5+R63.6+R71.8+R76.0+R76.8+R82.79+R91.8+S80.00XA+Y92.099+Z03.71+Z03.89+Z12.4+Z20.1+Z30.2+Z31.430+Z36.82+Z38.2+Z39.2+Z48.22+Z51.81+Z53.1+Z59.1+Z59.9+Z65.9+Z67.41+Z68.26+Z68.27+Z68.33+Z68.36+Z68.42+Z80.8+Z80.9+Z82.0+Z84.81+Z85.850+Z87.19+Z87.39+Z87.74+Z88.0+Z90.49+Z91.040+Z91.19+Z91.49+Z94.1,
                         rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_train, lambda=lambda[opt]))
  if(class(glm1)!="try-error"){  
    predictions[[i]] <- predict(glm_final, EMR_long_diags_test, type="response",s=lambda[opt])
    original[[i]]<-EMR_long_diags_test$Term
  }
  
}
save(predictions,original,file="predictions_diags.Rdata")

##After running
load("Data/predictions_diags.Rdata")
library(AUC)
auc<-NULL
for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}

tiff("AUC.tiff",res=300,w=2000,h=2000)
plot(roc(pred_train,factor(EMR_long_diags_train$Term)))
dev.off()


