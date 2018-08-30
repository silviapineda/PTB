rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB meds
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
load("Data/EMR_long_meds_multi.Rdata")

#demographics<-read.csv("Data/EMR_patients_patient_race.csv") 
##I am going to wait to Idit for the demo, meantime I will work only with the meds
# EMR_long_meds_multivariate$ID<-rownames(EMR_long_meds_multivariate)
# EMR_meds_demo<-merge(EMR_long_meds_multivariate,demographics,by="Patient_index")
# EMR_meds_demo<-subset(EMR_meds_demo,select = -c(ID,X,Term.y))
# colnames(EMR_meds_demo)[1]<-"Term"
# ##near Zero Variance correction
# id_nzv<-nearZeroVar(EMR_meds_demo,freqCut = 99/1,uniqueCut = 1)
# EMR_meds_demo_filter<-EMR_meds_demo[,-id_nzv] ##106

EMR_meds_demo<-EMR_long_meds_multivariate

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_meds_demo),collapse ="+")

########################
##### Downsampling ####
######################
set.seed(54)
for (i in 1:10){
  print(paste("down",i))
  ##Number of unique patient_index to obtain the downsampling
  EMR_long_meds_patient<-EMR_meds_demo[match(unique(EMR_meds_demo$Patient_index),
                                             EMR_meds_demo$Patient_index),c("Term","Patient_index")]
  
  dim_term<-table(EMR_long_meds_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_meds_patient[,"Term"])[2]
  id_down<-sample(c(1:dim_term),dim_ptb)
  EMR_down_term<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==0),][id_down,]
  EMR_down_ptb<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==1),]
  EMR_down<-rbind(EMR_down_term,EMR_down_ptb)
  id_down_data<-unlist(lapply(EMR_down$Patient_index, function(x) grep(x,EMR_long_meds_multivariate$Patient_index)))
  EMR_long_meds_multivariate_down<-EMR_long_meds_multivariate[id_down_data,]
  
  #########################
  ##  Run in the server  ##
  library(glmmLasso)
  
  ##Number of unique patient_index to obtain the train and test set
  EMR_long_meds_patient<-EMR_long_meds_multivariate_down[match(unique(EMR_long_meds_multivariate_down$Patient_index),
                                                               EMR_long_meds_multivariate_down$Patient_index),c("Term","Patient_index")]
  dim_term<-table(EMR_long_meds_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_meds_patient[,"Term"])[2]
  
  predictions<-list()
  original<-list()
  set.seed(54)
  for(j in 1:10){
    print(paste("Sample ", j,sep=""))
    id_test_data_PTB<-sample(c(1:dim_ptb),dim_ptb*0.2)
    EMR_test_PTB<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==1)[id_test_data_PTB],]
    id_test_data_Term<-sample(c(1:dim_term),dim_term*0.2)
    EMR_test_Term<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==0)[id_test_data_Term],]
    
    EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
    id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_meds_multivariate_down$Patient_index)))
    EMR_long_meds_test<-EMR_long_meds_multivariate_down[id_test_data,] 
    EMR_long_meds_train<-EMR_long_meds_multivariate_down[-(id_test_data),] 
    
    
    lambda <- seq(100,0,by=-5)
    family = binomial(link = logit)
    ################## First Simple Method ############################################
    ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
    
    BIC_vec<-rep(Inf,length(lambda))
    for(k in 1:length(lambda)){
      print(paste("Iteration ", k,sep=""))
      glm1 <- try(glmmLasso(Term ~ ACETAMINOPHEN..Acetaminophen.+ACYCLOVIR..Acyclovir.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+ASCORBIC.ACID..Ascorbic.Acid.+ASPIRIN..Aspirin.+AZITHROMYCIN..Azithromycin.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+Blood.Sugar..Blood.Glucose.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CLOTRIMAZOLE..Clotrimazole.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+DOCUSATE.SODIUM..Docusate.Sodium.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.SULFATE..ferrous.sulfate.+ferrous.sulfate.iron..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+GLUCOSE..Glucose.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+HYDROCORTISONE..Hydrocortisone.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORAZEPAM..Lorazepam.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+METFORMIN..Metformin.+METRONIDAZOLE..Metronidazole.+MICONAZOLE.NITRATE..Miconazole.Nitrate.+NIFEDIPINE..Nifedipine.+nitrofurantoin.macrocrystals..NITROFURANTOIN..MACROCRYSTALS.+ondansetron.hcl..Ondansetron.Hydrochloride.+PREDNISONE..Prednisone.+RANITIDINE..Ranitidine.+RHO.D.IMMUNE.GLOBULIN..Rho.D..Immune.Globulin.+Sennosides+SERTRALINE..Sertraline.+SIMETHICONE..Simethicone.+SODIUM.CHLORIDE..Sodium.Chloride.+Tetanus..tetanus.toxoid.vaccine..inactivated.+TRIAMCINOLONE.ACETONIDE..Triamcinolone.Acetonide.+Tuberculin.PPD..Purified.Protein.Derivative.of.Tuberculin.
                            ,rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_train, lambda=lambda[j]),silent=TRUE)  
      if(class(glm1)!="try-error"){  
        BIC_vec[j]<-glm1$bic
      }
    }
    
    opt<-which.min(BIC_vec)
    glm_final <- try(glmmLasso(Term ~ ACETAMINOPHEN..Acetaminophen.+ACYCLOVIR..Acyclovir.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+ASCORBIC.ACID..Ascorbic.Acid.+ASPIRIN..Aspirin.+AZITHROMYCIN..Azithromycin.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+Blood.Sugar..Blood.Glucose.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CLOTRIMAZOLE..Clotrimazole.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+DOCUSATE.SODIUM..Docusate.Sodium.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.SULFATE..ferrous.sulfate.+ferrous.sulfate.iron..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+GLUCOSE..Glucose.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+HYDROCORTISONE..Hydrocortisone.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORAZEPAM..Lorazepam.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+METFORMIN..Metformin.+METRONIDAZOLE..Metronidazole.+MICONAZOLE.NITRATE..Miconazole.Nitrate.+NIFEDIPINE..Nifedipine.+nitrofurantoin.macrocrystals..NITROFURANTOIN..MACROCRYSTALS.+ondansetron.hcl..Ondansetron.Hydrochloride.+PREDNISONE..Prednisone.+RANITIDINE..Ranitidine.+RHO.D.IMMUNE.GLOBULIN..Rho.D..Immune.Globulin.+Sennosides+SERTRALINE..Sertraline.+SIMETHICONE..Simethicone.+SODIUM.CHLORIDE..Sodium.Chloride.+Tetanus..tetanus.toxoid.vaccine..inactivated.+TRIAMCINOLONE.ACETONIDE..Triamcinolone.Acetonide.+Tuberculin.PPD..Purified.Protein.Derivative.of.Tuberculin.
                               ,rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_train, lambda=lambda[opt]))
    if(class(glm_final)!="try-error"){  
      summary(glm_final)
      
      predictions[[j]] <- predict(glm_final, EMR_long_meds_test, type="response",s=lambda[opt])
      original[[j]]<-EMR_long_meds_test$Term
    }
  }
  save(predictions,original,file=paste0("Results/predictions_meds_down_",i,".Rdata"))
}

####Run AUC
auc_mean<-NULL
for(i in 1:10){
  print(i)
  load(paste0("Results/predictions_meds_down_",i,".Rdata"))
  auc<-NULL
  for (j in 1:10){
    print(j)
    if(length(predictions[[j]])!=0){
      auc[j]<-auc(roc(predictions[[j]],factor(original[[j]])))
    }
  }
  auc_mean[i]<-mean(na.omit(auc))
}
