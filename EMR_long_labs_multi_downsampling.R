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

working_directory<-"/home/pinedasans/PTB/"
setwd(working_directory)

#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
load("Data/EMR_long_labs_multi_all.Rdata")

########################
##### Downsampling ####
######################
set.seed(54)
for (i in 1:10){
  print(paste("down",i))
  ##Number of unique patient_index to obtain the downsampling
  EMR_long_labs_patient<-EMR_long_labs_multivariate[match(unique(EMR_long_labs_multivariate$Patient_index),
                                                          EMR_long_labs_multivariate$Patient_index),c("Term","Patient_index")]
  
  dim_term<-table(EMR_long_labs_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_labs_patient[,"Term"])[2]
  id_down<-sample(c(1:dim_term),dim_ptb)
  EMR_down_term<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==0),][id_down,]
  EMR_down_ptb<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==1),]
  EMR_down<-rbind(EMR_down_term,EMR_down_ptb)
  id_down_data<-unlist(lapply(EMR_down$Patient_index, function(x) grep(x,EMR_long_labs_multivariate$Patient_index)))
  EMR_long_labs_multivariate_down<-EMR_long_labs_multivariate[id_down_data,]
  
  #########################
  ##  Run in the server  ##
  library(glmmLasso)
  
  ##Number of unique patient_index to obtain the train and test set
  EMR_long_labs_patient<-EMR_long_labs_multivariate_down[match(unique(EMR_long_labs_multivariate_down$Patient_index),
                                                               EMR_long_labs_multivariate_down$Patient_index),c("Term","Patient_index")]
  dim_term<-table(EMR_long_labs_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_labs_patient[,"Term"])[2]
  
  predictions<-list()
  original<-list()
  set.seed(54)
  for(j in 1:10){
    print(paste("Sample ", j,sep=""))
    id_test_data_PTB<-sample(c(1:dim_ptb),dim_ptb*0.2)
    EMR_test_PTB<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==1)[id_test_data_PTB],]
    id_test_data_Term<-sample(c(1:dim_term),dim_term*0.2)
    EMR_test_Term<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==0)[id_test_data_Term],]
    
    EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
    id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_labs_multivariate_down$Patient_index)))
    EMR_long_labs_test<-EMR_long_labs_multivariate_down[id_test_data,] 
    EMR_long_labs_train<-EMR_long_labs_multivariate_down[-(id_test_data),] 
    
    
    lambda <- seq(100,0,by=-5)
    family = binomial(link = logit)
    ################## First Simple Method ############################################
    ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
    
    BIC_vec<-rep(Inf,length(lambda))
    for(k in 1:length(lambda)){
      print(paste("Iteration ", k,sep=""))
      glm1 <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Banding.Resolution._order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Direct_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Biophysical.Profile.Score..of.8._order+C.Reactive.Protein_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Complement.C3..serum_order+Complement.C4..serum_order+Creat.per.Day..UR_order+Creatinine_order+
                              Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose...UA._order+Glucose..120.Mins_order+Glucose..180.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+HBV.Log.IU.mL_order+HBV.Real.Time_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A_order+
                              Hemoglobin.A1c_order+Hemoglobin.A2_order+Heparin.Level_order+Hours.Collected_order+IgG..serum_order+Imm.Gran..Left.Shift_order+Ketones..UA_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lipase_order+LVEDVi_order+LVEF.by.MOD.Bi.plane_order+LVESVi_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Metaphases.Analyzed._order+Metaphases.Counted._order+Metaphases.Karyotyped._order+Mins.Post.Gluc.Dose_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+Number.of.Cultures._order+OGT.Glucola.Dose_order+P.R.Interval_order+Parathormone_order+Parvo.Ab.B19.IgG_order+Parvo.Ab.B19.IgM_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+POCT.Glucose_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+
                              Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+RDW_order+Rubella.Antibody_order+RVVT.Seconds_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+T4..Total_order+Tacrolimus_order+Thyroid.Stimulating.Hormone_order+Thyroid.Stimulating.Immunoglobulin_order+Total.Volume.Collected_order+Triglycerides..serum_order+Urea.Nitrogen..Serum...Plasma_order+Uric.Acid..Serum...Plasma_order+Ventricular.Rate_order+Vitamin.B12_order+Vitamin.D..25.Hydroxy_order+WBC.Count_order+Weight.In.Kg_order,
                            rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[j]),silent=TRUE)  
      if(class(glm1)!="try-error"){  
        BIC_vec[j]<-glm1$bic
      }
    }
    
    opt<-which.min(BIC_vec)
    glm_final <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Banding.Resolution._order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Direct_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Biophysical.Profile.Score..of.8._order+C.Reactive.Protein_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Complement.C3..serum_order+Complement.C4..serum_order+Creat.per.Day..UR_order+Creatinine_order+
                             Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose...UA._order+Glucose..120.Mins_order+Glucose..180.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+HBV.Log.IU.mL_order+HBV.Real.Time_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A_order+
                             Hemoglobin.A1c_order+Hemoglobin.A2_order+Heparin.Level_order+Hours.Collected_order+IgG..serum_order+Imm.Gran..Left.Shift_order+Ketones..UA_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lipase_order+LVEDVi_order+LVEF.by.MOD.Bi.plane_order+LVESVi_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Metaphases.Analyzed._order+Metaphases.Counted._order+Metaphases.Karyotyped._order+Mins.Post.Gluc.Dose_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+Number.of.Cultures._order+OGT.Glucola.Dose_order+P.R.Interval_order+Parathormone_order+Parvo.Ab.B19.IgG_order+Parvo.Ab.B19.IgM_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+POCT.Glucose_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+
                             Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+RDW_order+Rubella.Antibody_order+RVVT.Seconds_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+T4..Total_order+Tacrolimus_order+Thyroid.Stimulating.Hormone_order+Thyroid.Stimulating.Immunoglobulin_order+Total.Volume.Collected_order+Triglycerides..serum_order+Urea.Nitrogen..Serum...Plasma_order+Uric.Acid..Serum...Plasma_order+Ventricular.Rate_order+Vitamin.B12_order+Vitamin.D..25.Hydroxy_order+WBC.Count_order+Weight.In.Kg_order,
                           rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[opt]))
    if(class(glm_final)!="try-error"){  
      summary(glm_final)
    
      predictions[[j]] <- predict(glm_final, EMR_long_labs_test, type="response",s=lambda[opt])
      original[[j]]<-EMR_long_labs_test$Term
    }
  }
  save(predictions,original,file=paste0("Results/predictions_labs_down_",i,".Rdata"))
}

