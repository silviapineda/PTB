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
### DESCRIP: Downsamplin 1:1 
###         Separate in training and testing 10 times and run glmmLasso
###         
###
### Author: Silvia Pineda
### Date: April, 2018
############################################################################################
library(lattice)
library(lme4)
library("RColorBrewer")
library(ggplot2)
library(caret)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)

#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
load("Data/EMR_long_diags_multi.Rdata")

#demographics<-read.csv("Data/EMR_patients_patient_race.csv") 
##I am going to wait to Idit for the demo, meantime I will work only with the diags
# EMR_long_diags_multivariate$ID<-rownames(EMR_long_diags_multivariate)
# EMR_diags_demo<-merge(EMR_long_diags_multivariate,demographics,by="Patient_index")
# EMR_diags_demo<-subset(EMR_diags_demo,select = -c(ID,X,Term.y))
# colnames(EMR_diags_demo)[1]<-"Term"
# ##near Zero Variance correction
# id_nzv<-nearZeroVar(EMR_diags_demo,freqCut = 99/1,uniqueCut = 1)
# EMR_diags_demo_filter<-EMR_diags_demo[,-id_nzv] ##106

EMR_diags_demo<-EMR_long_diags_multivariate

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_diags_demo),collapse ="+")

########################
##### Downsampling ####
######################
set.seed(54)
for (i in 1:10){
  print(paste("down",i))
  ##Number of unique patient_index to obtain the downsampling
  EMR_long_diags_patient<-EMR_diags_demo[match(unique(EMR_diags_demo$Patient_index),
                                             EMR_diags_demo$Patient_index),c("Term","Patient_index")]
  
  dim_term<-table(EMR_long_diags_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_diags_patient[,"Term"])[2]
  id_down<-sample(c(1:dim_term),dim_ptb)
  EMR_down_term<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==0),][id_down,]
  EMR_down_ptb<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==1),]
  EMR_down<-rbind(EMR_down_term,EMR_down_ptb)
  id_down_data<-unlist(lapply(EMR_down$Patient_index, function(x) grep(x,EMR_diags_demo$Patient_index)))
  EMR_long_diags_multivariate_down<-EMR_diags_demo[id_down_data,]
  
  #########################
  ##  Run in the server  ##
  library(glmmLasso)
  
  ##Number of unique patient_index to obtain the train and test set
  EMR_long_diags_patient<-EMR_long_diags_multivariate_down[match(unique(EMR_long_diags_multivariate_down$Patient_index),
                                                               EMR_long_diags_multivariate_down$Patient_index),c("Term","Patient_index")]
  dim_term<-table(EMR_long_diags_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_diags_patient[,"Term"])[2]
  
  predictions<-list()
  original<-list()
  set.seed(54)
  for(j in 1:10){
    print(paste("Sample ", j,sep=""))
    id_test_data_PTB<-sample(c(1:dim_ptb),dim_ptb*0.2)
    EMR_test_PTB<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==1)[id_test_data_PTB],]
    id_test_data_Term<-sample(c(1:dim_term),dim_term*0.2)
    EMR_test_Term<-EMR_long_diags_patient[which(EMR_long_diags_patient$Term==0)[id_test_data_Term],]
    
    EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
    id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_diags_multivariate_down$Patient_index)))
    EMR_long_diags_test<-EMR_long_diags_multivariate_down[id_test_data,] 
    EMR_long_diags_train<-EMR_long_diags_multivariate_down[-(id_test_data),] 
    
    
    lambda <- seq(100,0,by=-5)
    family = binomial(link = logit)
    ################## First Simple Method ############################################
    ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
    
    BIC_vec<-rep(Inf,length(lambda))
    for(k in 1:length(lambda)){
      print(paste("Iteration ", k,sep=""))
      glm1 <- try(glmmLasso(Term ~ D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.3+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O47.00+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891
                            ,rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_train, lambda=lambda[j]),silent=TRUE)  
      if(class(glm1)!="try-error"){  
        BIC_vec[j]<-glm1$bic
      }
    }
    
    opt<-which.min(BIC_vec)
    glm_final <- try(glmmLasso(Term ~ D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.3+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O47.00+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891
                               ,rnd = list(Patient_index=~1),family = family, data = EMR_long_diags_test, lambda=lambda[opt]))
    if(class(glm_final)!="try-error"){  
      summary(glm_final)
      
      predictions[[j]] <- predict(glm_final, EMR_long_diags_test, type="response",s=lambda[opt])
      original[[j]]<-EMR_long_diags_test$Term
    }
  }
  save(predictions,original,file=paste0("Results/predictions_diags_down_",i,".Rdata"))
}

####Run AUC
auc_mean<-NULL
for(i in 1:10){
  print(i)
  load(paste0("Results/predictions_diags_down_",i,".Rdata"))
  auc<-NULL
  for (j in 1:10){
    print(j)
    if(length(predictions[[j]])!=0){
      auc[j]<-auc(roc(predictions[[j]],factor(original[[j]])))
    }
  }
  auc_mean[i]<-mean(na.omit(auc))
}
