rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB diags
###
### CITATION: 
###
### PROCESS: Longitudinal data analysis. ENET for longitudinal data with all the data for the lab test
###           
### DESCRIP: Separate in training and testing 10 times and run glmmLasso
###         
###
### Author: Silvia Pineda
### Date: August, 2018
############################################################################################
library(lattice)
library(lme4)
library("RColorBrewer")
library(ggplot2)
library(caret)
library(glmmLasso)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)

#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
load("Data/EMR_long_diags_multi.Rdata")

demographics<-read.csv("Data/EMR_patients_patient_race_filtered.csv") 
demographics$Patient_Marital_Status<-gsub("RDP-Dissolved","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("RDP-LG SEP","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("Widowed","Unknown/Declined",demographics$Patient_Marital_Status)

demographics$Patient_Smoking_Status<-gsub("Heavy Tobacco Smoker","Current Every Day Smoker",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Never Assessed","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Unknown If Ever Smoked","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Smoker, Current Status Unknown","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)

EMR_long_diags_multivariate$ID<-rownames(EMR_long_diags_multivariate)
EMR_diags_demo<-merge(EMR_long_diags_multivariate,demographics,by="Patient_index")
EMR_diags_demo<-subset(EMR_diags_demo,select = -c(ID,X,Term.y))
colnames(EMR_diags_demo)[2]<-"Term"
# ##near Zero Variance correction
id_nzv<-nearZeroVar(EMR_diags_demo,freqCut = 99/1,uniqueCut = 1)
EMR_diags_demo_filter<-EMR_diags_demo[,-id_nzv] ##

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_diags_demo_filter),collapse ="+")

####################################
##  Run glmmLasso for the multi  ##
###################################

##Number of unique patient_index to obtain the train and test set
set.seed(54)
EMR_long_diags_patient<-EMR_diags_demo[match(unique(EMR_diags_demo_filter$Patient_index),EMR_diags_demo_filter$Patient_index)
                                     ,c("Term","Patient_index")]
table(EMR_long_diags_patient[,"Term"])
dim_term<-table(EMR_long_diags_patient[,"Term"])[1] ##8648 unique Term
dim_ptb<-table(EMR_long_diags_patient[,"Term"])[2]  ##414 unique PTB

predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:dim_ptb),0.2*dim_ptb) ##20% of PTB
  EMR_test_PTB<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:dim_term),0.2*dim_term) ##20% of term
  EMR_test_Term<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==0)[id_test_data_Term],]
  
  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_diags_demo_filter$Patient_index)))
  EMR_long_diags_test<-EMR_diags_demo_filter[id_test_data,] 
  EMR_long_diags_train<-EMR_diags_demo_filter[-(id_test_data),] 
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891
                          +Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                          ,rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891+Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                         ,rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_test, lambda=lambda[opt])
  summary(glm_final)
  
  predictions[[i]] <- predict(glm_final, EMR_long_diags_test, type="response",s=lambda[opt])
  original[[i]]<-EMR_long_diags_test$Term
}
save(predictions,original,file="Results/predictions_multi_diags.Rdata")



####
## Analysis Results
####
load("Results/predictions_multi_diags.Rdata")
library(AUC)
auc<-NULL

for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}

conf_matrix<-table(round(predictions[[i]]),factor(original[[i]]))
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)

tiff("Results/AUC_multi_diags.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC diags")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()
