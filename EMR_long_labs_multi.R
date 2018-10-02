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
library(glmmLasso)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)

#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
load("Data/EMR_long_labs_multi.Rdata")

demographics<-read.csv("Data/EMR_patients_patient_race_filtered.csv") 
demographics$Patient_Marital_Status<-gsub("RDP-Dissolved","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("RDP-LG SEP","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("Widowed","Unknown/Declined",demographics$Patient_Marital_Status)

demographics$Patient_Smoking_Status<-gsub("Heavy Tobacco Smoker","Current Every Day Smoker",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Never Assessed","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Unknown If Ever Smoked","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Smoker, Current Status Unknown","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)

EMR_long_labs_multivariate$ID<-rownames(EMR_long_labs_multivariate)
EMR_labs_demo<-merge(EMR_long_labs_multivariate,demographics,by="Patient_index")
EMR_labs_demo<-subset(EMR_labs_demo,select = -c(ID,X,Term.y))
colnames(EMR_labs_demo)[2]<-"Term"
##near Zero Variance correction
EMR_labs_demo_filter<-EMR_labs_demo

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_labs_demo_filter)[1:100],collapse ="+")
paste(colnames(EMR_labs_demo_filter)[101:200],collapse ="+")
paste(colnames(EMR_labs_demo_filter)[201:ncol(EMR_labs_demo_filter)],collapse ="+")
 

####################################
##  Run glmmLasso for the multi  ##
###################################

##Number of unique patient_index to obtain the train and test set
set.seed(54)
EMR_long_labs_patient<-EMR_labs_demo_filter[match(unique(EMR_labs_demo_filter$Patient_index),EMR_labs_demo_filter$Patient_index)
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
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_labs_demo_filter$Patient_index)))
  EMR_long_labs_test<-EMR_labs_demo_filter[id_test_data,] 
  EMR_long_labs_train<-EMR_labs_demo_filter[-(id_test_data),] 
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+X..Saturation_order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+AGE_order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Amylase..Serum...Plasma_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Banding.Resolution._order+Basophil.Abs.Count_order+Bicarbonate_order+Bile.Acids..Total_order+Bilirubin..Direct_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Biophysical.Profile.Score..of.8._order+BMI_order+C.Reactive.Protein_order+Calcium..Ionized..serum.plasma_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Cancer.Antigen.125_order+Carbon.Dioxide..Total_order+CD3.T.Cells.._order+CD3.T.Cells.Abs_order+CD4.CD8.Ratio_order+CD4.T.Cells.._order+CD4.T.Cells.Abs_order+CD8.T.Cells.._order+CD8.T.Cells.Abs_order+Chloride..Serum...Plasma_order+Colonies.Counted._order+Comments_order+Complement.C3..serum_order+Complement.C4..serum_order+Creat.per.Day..UR_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+Double.Stranded.DNA.Antibody_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Factor.VIII.Act..Average_order+FEF.25.75._order+FEF.25.75....PRED_order+FEF.25._order+FEF.25....PRED_order+FEF.50._order+FEF.50....PRED_order+FEF.75._order+FEF.75....PRED_order+Ferritin_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FEV1_order+FEV1...FVC_order+FEV1...FVC..PRED_order+FEV1..PRED_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Folate..RBC_order+Free.T3..Adult_order+Free.T4_order+FVC_order+FVC..PRED_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose...UA._order+Glucose..120.Mins_order+Glucose..180.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+HBV.Log.IU.mL_order+HBV.Real.Time_order+HEIGHT_order+Hematocrit_order+Hemoglobin_order+Hemoglobin..UA._order+Hemoglobin.A_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Hemoglobin.F_order+Heparin.Level_order+
                            Homocysteine..Total_order+Hours.Collected_order+ICD9.Code_order+Imm.Gran..Left.Shift_order+Iron..serum_order+Ketones..UA_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lactate..plasma_order+Lipase_order+LV.Mass.Index..BSA..by.M.Mode_order+LVEDVi_order+LVEF.by.MOD.Bi.plane_order+LVESVi_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Metaphases.Analyzed._order+Metaphases.Counted._order+Metaphases.Karyotyped._order+Methylmalonic.Acid..serum_order+Mins.Post.Gluc.Dose_order+Mitogen.Control_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+NIL.Control_order+Nonstress.Test_order+Number.of.Cultures._order+OGT.Glucola.Dose_order+Oxygen.Saturation_order+P.R.Interval_order+Parathormone_order+Parvo.Ab.B19.IgG_order+Parvo.Ab.B19.IgM_order+PCO2_order+PEAK.FLOW_order+PEAK.FLOW..PRED_order+pH..Blood_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+Platelets...Units.Ready_order+PO2_order+POCT.Glucose_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Potassium..whole.blood_order+Prealbumin_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+RDW_order+Retic.Count..Flow.Cytometry_order+Ristocetin.Cofactor_order+Rubella.Antibody_order+RVVT.Seconds_order+Sedimentation.Rate_order+Sedimentation.Rate..MB._order+Sodium..Serum...Plasma_order+Sodium..whole.blood_order+Specific.Gravity_order+T3..Total_order+T4..Total_order+Tacrolimus_order+TB.Antigen_order+Thyroid.Stimulating.Hormone_order+Thyroid.Stimulating.Immunoglobulin_order+Total.Volume.Collected_order+Transferrin_order+Triglycerides..serum_order+Urea.Nitrogen..Serum...Plasma_order+Uric.Acid..Serum...Plasma_order+Urobilinogen_order+Vancomycin_order+Ventricular.Rate_order+Vitamin.B12_order+Vitamin.D..25.Hydroxy_order+von.Willebrand.Factor.Antigen_order+WBC.Count_order+WEIGHT_order+Weight.In.Kg_order+X.Int.l.Normaliz.Ratio.+X..Saturation+Activated.Partial.Thromboplastin.Time+AFI..cm.+AGE+Alanine.transaminase+Albumin..Random+Albumin..Serum...Plasma+Alkaline.Phosphatase+
                            Amylase..Serum...Plasma+Anion.Gap+Aspartate.transaminase+Atrial.Rate+Banding.Resolution.+Basophil.Abs.Count+Bicarbonate+Bile.Acids..Total+Bilirubin..Direct+Bilirubin..Total+Biophysical.Profile.Score..of.10.+Biophysical.Profile.Score..of.8.+BMI+C.Reactive.Protein+Calcium..Ionized..serum.plasma+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Cancer.Antigen.125+Carbon.Dioxide..Total+CD3.T.Cells..+CD3.T.Cells.Abs+CD4.CD8.Ratio+CD4.T.Cells..+CD4.T.Cells.Abs+CD8.T.Cells..+CD8.T.Cells.Abs+Chloride..Serum...Plasma+Colonies.Counted.+Comments+Complement.C3..serum+Complement.C4..serum+Creat.per.Day..UR+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+Double.Stranded.DNA.Antibody+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Factor.VIII.Act..Average+FEF.25.75.+FEF.25.75....PRED+FEF.25.+FEF.25....PRED+FEF.50.+FEF.50....PRED+FEF.75.+FEF.75....PRED+Ferritin+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FEV1+FEV1...FVC+FEV1...FVC..PRED+FEV1..PRED+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Folate..RBC+Free.T3..Adult+Free.T4+FVC+FVC..PRED+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose.Loading.Screen+Glucose...UA.+Glucose..120.Mins+Glucose..180.Mins+Glucose..60.Mins+Glucose..fasting+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..Fasting..Pregnant..3.hr.GTT+Glucose..meter.download+Glucose..non.fasting+HBV.Log.IU.mL+HBV.Real.Time+HEIGHT+Hematocrit+Hemoglobin+Hemoglobin..UA.+Hemoglobin.A+Hemoglobin.A1c+Hemoglobin.A2+Hemoglobin.F+Heparin.Level+Homocysteine..Total+Hours.Collected+ICD9.Code+Imm.Gran..Left.Shift+Iron..serum+Ketones..UA+Lactate.Dehydrogenase..Serum...Plasma+Lactate..plasma+Lipase+LV.Mass.Index..BSA..by.M.Mode+LVEDVi+LVEF.by.MOD.Bi.plane+LVESVi+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+MCH+MCHC+MCV+Metaphases.Analyzed.+Metaphases.Counted.+Metaphases.Karyotyped.+Methylmalonic.Acid..serum+Mins.Post.Gluc.Dose+Mitogen.Control+Monocyte.Abs.Count+Neutrophil.Absolute.Count+NIL.Control+Nonstress.Test+Number.of.Cultures.+OGT.Glucola.Dose+Oxygen.Saturation+P.R.Interval+Parathormone+Parvo.Ab.B19.IgG+Parvo.Ab.B19.IgM+PCO2+PEAK.FLOW+PEAK.FLOW..PRED+pH..Blood+pH..UA+Phosphorus..Serum...Plasma+Platelet.Count+Platelets...Units.Ready+PO2+POCT.Glucose+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Potassium..whole.blood+Prealbumin+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+Protein..Total..Serum...Plasma+Protein..UA+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+RBCs...Units.Ready+RDW+Retic.Count..Flow.Cytometry+Ristocetin.Cofactor+Rubella.Antibody+RVVT.Seconds+Sedimentation.Rate+Sedimentation.Rate..MB.+Sodium..Serum...Plasma+Sodium..whole.blood+Specific.Gravity+T3..Total+T4..Total+Tacrolimus+TB.Antigen+Thyroid.Stimulating.Hormone+Thyroid.Stimulating.Immunoglobulin+Total.Volume.Collected+Transferrin+Triglycerides..serum+Urea.Nitrogen..Serum...Plasma+Uric.Acid..Serum...Plasma+Urobilinogen+Vancomycin+Ventricular.Rate+Vitamin.B12+Vitamin.D..25.Hydroxy+von.Willebrand.Factor.Antigen+WBC.Count+
                            WEIGHT+Weight.In.Kg+Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                          ,rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+X..Saturation_order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+AGE_order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Amylase..Serum...Plasma_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Banding.Resolution._order+Basophil.Abs.Count_order+Bicarbonate_order+Bile.Acids..Total_order+Bilirubin..Direct_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Biophysical.Profile.Score..of.8._order+BMI_order+C.Reactive.Protein_order+Calcium..Ionized..serum.plasma_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Cancer.Antigen.125_order+Carbon.Dioxide..Total_order+CD3.T.Cells.._order+CD3.T.Cells.Abs_order+CD4.CD8.Ratio_order+CD4.T.Cells.._order+CD4.T.Cells.Abs_order+CD8.T.Cells.._order+CD8.T.Cells.Abs_order+Chloride..Serum...Plasma_order+Colonies.Counted._order+Comments_order+Complement.C3..serum_order+Complement.C4..serum_order+Creat.per.Day..UR_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+Double.Stranded.DNA.Antibody_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Factor.VIII.Act..Average_order+FEF.25.75._order+FEF.25.75....PRED_order+FEF.25._order+FEF.25....PRED_order+FEF.50._order+FEF.50....PRED_order+FEF.75._order+FEF.75....PRED_order+Ferritin_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FEV1_order+FEV1...FVC_order+FEV1...FVC..PRED_order+FEV1..PRED_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Folate..RBC_order+Free.T3..Adult_order+Free.T4_order+FVC_order+FVC..PRED_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose.Loading.Screen_order+Glucose...UA._order+Glucose..120.Mins_order+Glucose..180.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+HBV.Log.IU.mL_order+HBV.Real.Time_order+HEIGHT_order+Hematocrit_order+Hemoglobin_order+Hemoglobin..UA._order+Hemoglobin.A_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Hemoglobin.F_order+Heparin.Level_order+
                           Homocysteine..Total_order+Hours.Collected_order+ICD9.Code_order+Imm.Gran..Left.Shift_order+Iron..serum_order+Ketones..UA_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lactate..plasma_order+Lipase_order+LV.Mass.Index..BSA..by.M.Mode_order+LVEDVi_order+LVEF.by.MOD.Bi.plane_order+LVESVi_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Metaphases.Analyzed._order+Metaphases.Counted._order+Metaphases.Karyotyped._order+Methylmalonic.Acid..serum_order+Mins.Post.Gluc.Dose_order+Mitogen.Control_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+NIL.Control_order+Nonstress.Test_order+Number.of.Cultures._order+OGT.Glucola.Dose_order+Oxygen.Saturation_order+P.R.Interval_order+Parathormone_order+Parvo.Ab.B19.IgG_order+Parvo.Ab.B19.IgM_order+PCO2_order+PEAK.FLOW_order+PEAK.FLOW..PRED_order+pH..Blood_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+Platelets...Units.Ready_order+PO2_order+POCT.Glucose_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Potassium..whole.blood_order+Prealbumin_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+RDW_order+Retic.Count..Flow.Cytometry_order+Ristocetin.Cofactor_order+Rubella.Antibody_order+RVVT.Seconds_order+Sedimentation.Rate_order+Sedimentation.Rate..MB._order+Sodium..Serum...Plasma_order+Sodium..whole.blood_order+Specific.Gravity_order+T3..Total_order+T4..Total_order+Tacrolimus_order+TB.Antigen_order+Thyroid.Stimulating.Hormone_order+Thyroid.Stimulating.Immunoglobulin_order+Total.Volume.Collected_order+Transferrin_order+Triglycerides..serum_order+Urea.Nitrogen..Serum...Plasma_order+Uric.Acid..Serum...Plasma_order+Urobilinogen_order+Vancomycin_order+Ventricular.Rate_order+Vitamin.B12_order+Vitamin.D..25.Hydroxy_order+von.Willebrand.Factor.Antigen_order+WBC.Count_order+WEIGHT_order+Weight.In.Kg_order+X.Int.l.Normaliz.Ratio.+X..Saturation+Activated.Partial.Thromboplastin.Time+AFI..cm.+AGE+Alanine.transaminase+Albumin..Random+Albumin..Serum...Plasma+Alkaline.Phosphatase+
                           Amylase..Serum...Plasma+Anion.Gap+Aspartate.transaminase+Atrial.Rate+Banding.Resolution.+Basophil.Abs.Count+Bicarbonate+Bile.Acids..Total+Bilirubin..Direct+Bilirubin..Total+Biophysical.Profile.Score..of.10.+Biophysical.Profile.Score..of.8.+BMI+C.Reactive.Protein+Calcium..Ionized..serum.plasma+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Cancer.Antigen.125+Carbon.Dioxide..Total+CD3.T.Cells..+CD3.T.Cells.Abs+CD4.CD8.Ratio+CD4.T.Cells..+CD4.T.Cells.Abs+CD8.T.Cells..+CD8.T.Cells.Abs+Chloride..Serum...Plasma+Colonies.Counted.+Comments+Complement.C3..serum+Complement.C4..serum+Creat.per.Day..UR+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+Double.Stranded.DNA.Antibody+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Factor.VIII.Act..Average+FEF.25.75.+FEF.25.75....PRED+FEF.25.+FEF.25....PRED+FEF.50.+FEF.50....PRED+FEF.75.+FEF.75....PRED+Ferritin+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FEV1+FEV1...FVC+FEV1...FVC..PRED+FEV1..PRED+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Folate..RBC+Free.T3..Adult+Free.T4+FVC+FVC..PRED+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose.Loading.Screen+Glucose...UA.+Glucose..120.Mins+Glucose..180.Mins+Glucose..60.Mins+Glucose..fasting+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..Fasting..Pregnant..3.hr.GTT+Glucose..meter.download+Glucose..non.fasting+HBV.Log.IU.mL+HBV.Real.Time+HEIGHT+Hematocrit+Hemoglobin+Hemoglobin..UA.+Hemoglobin.A+Hemoglobin.A1c+Hemoglobin.A2+Hemoglobin.F+Heparin.Level+Homocysteine..Total+Hours.Collected+ICD9.Code+Imm.Gran..Left.Shift+Iron..serum+Ketones..UA+Lactate.Dehydrogenase..Serum...Plasma+Lactate..plasma+Lipase+LV.Mass.Index..BSA..by.M.Mode+LVEDVi+LVEF.by.MOD.Bi.plane+LVESVi+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+MCH+MCHC+MCV+Metaphases.Analyzed.+Metaphases.Counted.+Metaphases.Karyotyped.+Methylmalonic.Acid..serum+Mins.Post.Gluc.Dose+Mitogen.Control+Monocyte.Abs.Count+Neutrophil.Absolute.Count+NIL.Control+Nonstress.Test+Number.of.Cultures.+OGT.Glucola.Dose+Oxygen.Saturation+P.R.Interval+Parathormone+Parvo.Ab.B19.IgG+Parvo.Ab.B19.IgM+PCO2+PEAK.FLOW+PEAK.FLOW..PRED+pH..Blood+pH..UA+Phosphorus..Serum...Plasma+Platelet.Count+Platelets...Units.Ready+PO2+POCT.Glucose+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Potassium..whole.blood+Prealbumin+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+Protein..Total..Serum...Plasma+Protein..UA+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+RBCs...Units.Ready+RDW+Retic.Count..Flow.Cytometry+Ristocetin.Cofactor+Rubella.Antibody+RVVT.Seconds+Sedimentation.Rate+Sedimentation.Rate..MB.+Sodium..Serum...Plasma+Sodium..whole.blood+Specific.Gravity+T3..Total+T4..Total+Tacrolimus+TB.Antigen+Thyroid.Stimulating.Hormone+Thyroid.Stimulating.Immunoglobulin+Total.Volume.Collected+Transferrin+Triglycerides..serum+Urea.Nitrogen..Serum...Plasma+Uric.Acid..Serum...Plasma+Urobilinogen+Vancomycin+Ventricular.Rate+Vitamin.B12+Vitamin.D..25.Hydroxy+von.Willebrand.Factor.Antigen+WBC.Count+
                           WEIGHT+Weight.In.Kg+Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                         ,rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_test, lambda=lambda[opt])
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
