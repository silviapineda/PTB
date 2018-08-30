rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB labs
###
### CITATION: 
###
### PROCESS: Longitudinal data analysis. ENET for longitudinal data with all the data for the lab test
###           
### DESCRIP: Separate in training and testing 10 times and run glmmLasso
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
load("Data/EMR_long_labs_multi.Rdata")

#demographics<-read.csv("Data/EMR_patients_patient_race.csv") 
##I am going to wait to Idit for the demo, meantime I will work only with the labs
EMR_long_labs_multivariate$ID<-rownames(EMR_long_labs_multivariate)
EMR_labs_demo<-merge(EMR_long_labs_multivariate,demographics,by="Patient_index")
EMR_labs_demo<-subset(EMR_labs_demo,select = -c(ID,X,Term.y))
colnames(EMR_labs_demo)[1]<-"Term"
##near Zero Variance correction
id_nzv<-nearZeroVar(EMR_labs_demo,freqCut = 99/1,uniqueCut = 1)
EMR_labs_demo_filter<-EMR_labs_demo[,-id_nzv] ##106

EMR_labs_demo<-EMR_long_labs_multivariate

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_labs_demo)[1:30],collapse ="+")
paste(colnames(EMR_labs_demo)[31:60],collapse ="+")
paste(colnames(EMR_labs_demo)[61:ncol(EMR_labs_demo)],collapse ="+")
 

####################################
##  Run glmmLasso for the multi  ##
###################################

library(glmmLasso)

##Number of unique patient_index to obtain the train and test set
set.seed(54)
EMR_long_labs_patient<-EMR_labs_demo[match(unique(EMR_labs_demo$Patient_index),EMR_labs_demo$Patient_index)
                                                  ,c("Term","Patient_index")]
table(EMR_long_labs_patient[,"Term"])
dim_term<-table(EMR_long_labs_patient[,"Term"])[1] ##5005 unique Term
dim_ptb<-table(EMR_long_labs_patient[,"Term"])[2]  ##324 unique PTB

predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:dim_ptb),0.2*dim_ptb) ##20% of PTB
  EMR_test_PTB<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:dim_term),0.2*dim_term) ##20% of term
  EMR_test_Term<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==0)[id_test_data_Term],]
  
  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_labs_demo$Patient_index)))
  EMR_long_labs_test<-EMR_labs_demo[id_test_data,] 
  EMR_long_labs_train<-EMR_labs_demo[-(id_test_data),] 
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+
                            Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose..120.Mins_order+Glucose..60.Mins_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A1c_order+Hours.Collected_order+Lymphocyte.Abs.Cnt_order+MCH_order+MCHC_order+MCV_order+Mins.Post.Gluc.Dose_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+OGT.Glucola.Dose_order+P.R.Interval_order+pH..UA_order+Platelet.Count_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+
                            Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Urea.Nitrogen..Serum...Plasma_order+Ventricular.Rate_order+WBC.Count_order+X.Int.l.Normaliz.Ratio.+Activated.Partial.Thromboplastin.Time+AFI..cm.+Alanine.transaminase+Albumin..Serum...Plasma+Alkaline.Phosphatase+Anion.Gap+Aspartate.transaminase+Atrial.Rate+Basophil.Abs.Count+Bile.Acids..Total+Bilirubin..Total+Biophysical.Profile.Score..of.10.+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Carbon.Dioxide..Total+Chloride..Serum...Plasma+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Free.T4+Glu.Tol.Post.Glucola+Glucose.Loading.Screen+Glucose..120.Mins+Glucose..60.Mins+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..Fasting..Pregnant..3.hr.GTT+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin+Hemoglobin.A1c+Hours.Collected+Lymphocyte.Abs.Cnt+MCH+MCHC+MCV+Mins.Post.Gluc.Dose+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Nonstress.Test+OGT.Glucola.Dose+P.R.Interval+pH..UA+Platelet.Count+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Urea.Nitrogen..Serum...Plasma+Ventricular.Rate+WBC.Count
                          ,rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+
                           Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose..120.Mins_order+Glucose..60.Mins_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A1c_order+Hours.Collected_order+Lymphocyte.Abs.Cnt_order+MCH_order+MCHC_order+MCV_order+Mins.Post.Gluc.Dose_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+OGT.Glucola.Dose_order+P.R.Interval_order+pH..UA_order+Platelet.Count_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+
                           Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Urea.Nitrogen..Serum...Plasma_order+Ventricular.Rate_order+WBC.Count_order+X.Int.l.Normaliz.Ratio.+Activated.Partial.Thromboplastin.Time+AFI..cm.+Alanine.transaminase+Albumin..Serum...Plasma+Alkaline.Phosphatase+Anion.Gap+Aspartate.transaminase+Atrial.Rate+Basophil.Abs.Count+Bile.Acids..Total+Bilirubin..Total+Biophysical.Profile.Score..of.10.+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Carbon.Dioxide..Total+Chloride..Serum...Plasma+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Free.T4+Glu.Tol.Post.Glucola+Glucose.Loading.Screen+Glucose..120.Mins+Glucose..60.Mins+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..Fasting..Pregnant..3.hr.GTT+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin+Hemoglobin.A1c+Hours.Collected+Lymphocyte.Abs.Cnt+MCH+MCHC+MCV+Mins.Post.Gluc.Dose+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Nonstress.Test+OGT.Glucola.Dose+P.R.Interval+pH..UA+Platelet.Count+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Urea.Nitrogen..Serum...Plasma+Ventricular.Rate+WBC.Count
                         ,rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[opt])
  summary(glm_final)
  
  predictions[[i]] <- predict(glm_final, EMR_long_labs_test, type="response",s=lambda[opt])
  original[[i]]<-EMR_long_labs_test$Term
}
save(predictions,original,file="Results/predictions_multi_labs.Rdata")



####
## Analysis Results
####
load("Results/predictions_multi_labs.Rdata")
library(AUC)
auc<-NULL

for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}

conf_matrix<-table(round(predictions[[i]]),factor(original[[i]]))
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)

tiff("Results/AUC_multi_labs.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC labs")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()
