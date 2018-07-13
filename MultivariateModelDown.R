rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB Diags
###
### CITATION: 
###
### PROCESS: Adjust multivariable model with MEDS+LABS+DIAGS+Demographics (downsampling)
###           
### DESCRIP: Analysis of EMR data for PTB
###         
###
### Author: Silvia Pineda
### Date: July, 2018
############################################################################################
library("RColorBrewer")
library(ggplot2)
library("glmmLasso")
working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)


#########################
##  Run in the server  ##
library(glmmLasso)
load("Data/EMR_multi_all.Rdata")

# ##To extract the names of the variables
paste(colnames(EMR_labs_diags_meds_demo_filter)[1:100],collapse ="+")
paste(colnames(EMR_labs_diags_meds_demo_filter)[101:200],collapse ="+")
paste(colnames(EMR_labs_diags_meds_demo_filter)[201:300],collapse ="+")
paste(colnames(EMR_labs_diags_meds_demo_filter)[301:400],collapse ="+")
paste(colnames(EMR_labs_diags_meds_demo_filter)[401:ncol(EMR_labs_diags_meds_demo_filter)],collapse ="+")

########################
##### Downsampling ####
######################
set.seed(54)
for (i in 1:10){
  print(paste("down",i))
  ##Number of unique patient_index to obtain the downsampling
  EMR_long_patient<-EMR_labs_diags_meds_demo_filter[match(unique(EMR_labs_diags_meds_demo_filter$Patient_index),
                                                          EMR_labs_diags_meds_demo_filter$Patient_index),c("Term","Patient_index")]
  
  dim_term<-table(EMR_long_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_patient[,"Term"])[2]
  id_down<-sample(c(1:dim_term),dim_ptb*2)
  EMR_down_term<-EMR_long_patient[which(EMR_long_patient$Term==0),][id_down,]
  EMR_down_ptb<-EMR_long_patient[which(EMR_long_patient$Term==1),]
  EMR_down<-rbind(EMR_down_term,EMR_down_ptb)
  id_down_data<-unlist(lapply(EMR_down$Patient_index, function(x) grep(x,EMR_labs_diags_meds_demo_filter$Patient_index)))
  EMR_long_multivariate_down<-EMR_labs_diags_meds_demo_filter[id_down_data,] 
  
  ##Number of unique patient_index to obtain the train and test set
  EMR_long_patient<-EMR_long_multivariate_down[match(unique(EMR_long_multivariate_down$Patient_index),
                                                          EMR_long_multivariate_down$Patient_index),c("Term","Patient_index")]
  dim_term<-table(EMR_long_patient[,"Term"])[1]
  dim_ptb<-table(EMR_long_patient[,"Term"])[2]
  
  predictions<-list()
  original<-list()
  set.seed(54)
  
  for(j in 1:10){
    print(paste("Sample ", j,sep=""))
    id_test_data_PTB<-sample(c(1:dim_ptb),dim_ptb*0.2)
    EMR_test_PTB<-EMR_long_patient[which(EMR_long_patient$Term==1)[id_test_data_PTB],]
    id_test_data_Term<-sample(c(1:dim_term),dim_term*0.2)
    EMR_test_Term<-EMR_long_patient[which(EMR_long_patient$Term==0)[id_test_data_Term],]
    
    EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
    id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_multivariate_down$Patient_index)))
    EMR_long_test<-EMR_long_multivariate_down[id_test_data,] 
    EMR_long_train<-EMR_long_multivariate_down[-(id_test_data),] 
    
    # tab<-NULL
    # for (i in 1:ncol(EMR_long_train)){
    #   tab[i]<-dim(table(EMR_long_train[,i]))
    # }
    # 
    lambda <- seq(100,0,by=-5)
    family = binomial(link = logit)
    ################## First Simple Method ############################################
    ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
    BIC_vec<-rep(Inf,length(lambda))
  
    for(k in 1:length(lambda)){
      print(paste("Iteration ", k,sep=""))
      glm1 <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Complement.C3..serum_order+Creat.per.Day..UR_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Gamma.Glutamyl.Transpeptidase_order+Glucose.Loading.Screen_order+Glucose..120.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Hours.Collected_order+Imm.Gran..Left.Shift_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+OGT.Glucola.Dose_order+P.R.Interval_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Tacrolimus_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Urea.Nitrogen..Serum...Plasma_order+Ventricular.Rate_order+Vitamin.D..25.Hydroxy_order+WBC.Count_order+Weight.In.Kg_order+X.Int.l.Normaliz.Ratio.+Activated.Partial.Thromboplastin.Time+AFI..cm.+Alanine.transaminase+Albumin..Random+Albumin..Serum...Plasma+Alkaline.Phosphatase+
                              Anion.Gap+Aspartate.transaminase+Atrial.Rate+Basophil.Abs.Count+Bile.Acids..Total+Bilirubin..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Carbon.Dioxide..Total+Chloride..Serum...Plasma+Complement.C3..serum+Creat.per.Day..UR+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Free.T4+Gamma.Glutamyl.Transpeptidase+Glucose.Loading.Screen+Glucose..120.Mins+Glucose..60.Mins+Glucose..fasting+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin+Hemoglobin.A+Hemoglobin.A1c+Hemoglobin.A2+Hours.Collected+Imm.Gran..Left.Shift+Lactate.Dehydrogenase..Serum...Plasma+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+MCH+MCHC+MCV+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Nonstress.Test+OGT.Glucola.Dose+P.R.Interval+pH..UA+Phosphorus..Serum...Plasma+Platelet.Count+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+Protein..Total..Serum...Plasma+Protein..UA+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+RBCs...Units.Ready+Sodium..Serum...Plasma+Specific.Gravity+Tacrolimus+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Urea.Nitrogen..Serum...Plasma+Ventricular.Rate+Vitamin.D..25.Hydroxy+WBC.Count+Weight.In.Kg+B00.9+B37.3+D25.9+D63.1+D64.9+E03.9+E07.9+E11.65+E11.9+E66.9+F19.10+F32.9+F41.1+F41.8+F41.9+G43.909+I10+
                              I25.10+I42.9+I49.8+I82.409+IMO0001+IMO0002+J45.909+K21.9+K50.90+K83.1+L93.0+M32.9+M54.9+M89.9+N12+N18.9+N39.0+N76.0+N93.9+O09.212+O09.213+O09.219+O09.299+O09.513+O09.519+O09.522+O09.523+O09.529+O09.819+O09.892+O09.893+O09.899+O09.90+O09.92+O09.93+O10.019+O13.3+O13.9+O16.9+O22.30+O23.90+O24.410+O24.414+O24.419+O24.919+O26.613+O26.619+O26.843+O26.849+O26.853+O26.859+O26.873+O26.879+O26.893+O28.9+O32.1XX0+O34.21+O34.219+O34.30+O34.40+O35.0XX0+O35.8XX0+O35.9XX0+O35.9XX1+O36.0130+O36.0131+O36.0990+O36.5930+O36.5990+O36.8130+O36.8190+O36.8390+O36.8990+O40.3XX0+O40.9XX0+O41.00X0+O43.899+O43.93+O44.00+O44.03+O44.10+O45.90+O46.8X9+O46.90+O46.93+O47.00+O47.03+O47.9+O62.9+O76+O98.513+O98.519+O99.012+O99.013+O99.019+O99.119+O99.210+O99.213+O99.280+O99.282+
                              O99.283+O99.340+O99.343+O99.413+O99.419+O99.513+O99.613+O99.810+O99.820+O99.89+Q24.9+Q28.9+Q79.1+R00.0+R03.0+R10.9+R19.7+R21+R51+R76.11+R76.8+R80.9+R82.71+Z22.330+Z23+Z29.8+Z33.1+Z34.00+Z34.02+Z34.03+Z34.80+Z34.83+Z34.90+Z34.93+Z36+Z36.3+Z41.8+Z67.91+Z78.9+Z79.4+Z82.49+Z83.3+Z86.19+Z87.440+Z87.891+Z88.0+Z96.41+Z98.890+Z98.891+ACETAMINOPHEN..Acetaminophen.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+AMPICILLIN..Ampicillin.+ASPIRIN..Aspirin.+BECLOMETHASONE.DIPROPIONATE..Beclomethasone.Dipropionate.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+BUPIVACAINE..Bupivacaine.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CEFAZOLIN..Cefazolin.+CEFTRIAXONE..Ceftriaxone.+CEPHALEXIN..Cephalexin.+CETIRIZINE..Cetirizine.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+DOCUSATE.SODIUM..Docusate.Sodium.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+EPINEPHRINE..Epinephrine.+ERYTHROMYCIN..Erythromycin.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.SULFATE..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+FOLIC.ACID..Folic.Acid.+GLUCOSE..Glucose.+Glycol..polyethylene..Polyethylene.Glycols.+HEPARIN..Heparin.+HYDRALAZINE..Hydralazine.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+hydrocortisone.sod.succinate..Hydrocortisone.sodium.succinate.+HYDROMORPHONE..Hydromorphone.+HYDROXYPROGESTERONE.CAPROATE..hydroxyprogesterone.caproate..USP..+INDOMETHACIN..Indomethacin.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LANSOPRAZOLE..lansoprazole.+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORAZEPAM..Lorazepam.+MAGNESIUM.OXIDE..Magnesium.Oxide.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Makena+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+METFORMIN..Metformin.+METOCLOPRAMIDE..Metoclopramide.+MIDAZOLAM..Midazolam.+MORPHINE..Morphine.
                            ,rnd = list(Patient_index=~1),family = family, data = EMR_long_train, lambda=lambda[k]),silent=TRUE)
      if(class(glm1)!="try-error"){  
        BIC_vec[j]<-glm1$bic
      }
    }
  
    opt<-which.min(BIC_vec)
    glm_final <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Random_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Complement.C3..serum_order+Creat.per.Day..UR_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Gamma.Glutamyl.Transpeptidase_order+Glucose.Loading.Screen_order+Glucose..120.Mins_order+Glucose..60.Mins_order+Glucose..fasting_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Hours.Collected_order+Imm.Gran..Left.Shift_order+Lactate.Dehydrogenase..Serum...Plasma_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+MCH_order+MCHC_order+MCV_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+OGT.Glucola.Dose_order+P.R.Interval_order+pH..UA_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+Protein..Total..Serum...Plasma_order+Protein..UA_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+RBCs...Units.Ready_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Tacrolimus_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Urea.Nitrogen..Serum...Plasma_order+Ventricular.Rate_order+Vitamin.D..25.Hydroxy_order+WBC.Count_order+Weight.In.Kg_order+X.Int.l.Normaliz.Ratio.+Activated.Partial.Thromboplastin.Time+AFI..cm.+Alanine.transaminase+Albumin..Random+Albumin..Serum...Plasma+Alkaline.Phosphatase+
                                 Anion.Gap+Aspartate.transaminase+Atrial.Rate+Basophil.Abs.Count+Bile.Acids..Total+Bilirubin..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Carbon.Dioxide..Total+Chloride..Serum...Plasma+Complement.C3..serum+Creat.per.Day..UR+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FHR.Baseline+Fibrinogen..Functional+Fluid.Volume+Free.T4+Gamma.Glutamyl.Transpeptidase+Glucose.Loading.Screen+Glucose..120.Mins+Glucose..60.Mins+Glucose..fasting+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin+Hemoglobin.A+Hemoglobin.A1c+Hemoglobin.A2+Hours.Collected+Imm.Gran..Left.Shift+Lactate.Dehydrogenase..Serum...Plasma+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+MCH+MCHC+MCV+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Nonstress.Test+OGT.Glucola.Dose+P.R.Interval+pH..UA+Phosphorus..Serum...Plasma+Platelet.Count+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+Protein..Total..Serum...Plasma+Protein..UA+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+RBCs...Units.Ready+Sodium..Serum...Plasma+Specific.Gravity+Tacrolimus+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Urea.Nitrogen..Serum...Plasma+Ventricular.Rate+Vitamin.D..25.Hydroxy+WBC.Count+Weight.In.Kg+B00.9+B37.3+D25.9+D63.1+D64.9+E03.9+E07.9+E11.65+E11.9+E66.9+F19.10+F32.9+F41.1+F41.8+F41.9+G43.909+I10+
                                 I25.10+I42.9+I49.8+I82.409+IMO0001+IMO0002+J45.909+K21.9+K50.90+K83.1+L93.0+M32.9+M54.9+M89.9+N12+N18.9+N39.0+N76.0+N93.9+O09.212+O09.213+O09.219+O09.299+O09.513+O09.519+O09.522+O09.523+O09.529+O09.819+O09.892+O09.893+O09.899+O09.90+O09.92+O09.93+O10.019+O13.3+O13.9+O16.9+O22.30+O23.90+O24.410+O24.414+O24.419+O24.919+O26.613+O26.619+O26.843+O26.849+O26.853+O26.859+O26.873+O26.879+O26.893+O28.9+O32.1XX0+O34.21+O34.219+O34.30+O34.40+O35.0XX0+O35.8XX0+O35.9XX0+O35.9XX1+O36.0130+O36.0131+O36.0990+O36.5930+O36.5990+O36.8130+O36.8190+O36.8390+O36.8990+O40.3XX0+O40.9XX0+O41.00X0+O43.899+O43.93+O44.00+O44.03+O44.10+O45.90+O46.8X9+O46.90+O46.93+O47.00+O47.03+O47.9+O62.9+O76+O98.513+O98.519+O99.012+O99.013+O99.019+O99.119+O99.210+O99.213+O99.280+O99.282+
                                 O99.283+O99.340+O99.343+O99.413+O99.419+O99.513+O99.613+O99.810+O99.820+O99.89+Q24.9+Q28.9+Q79.1+R00.0+R03.0+R10.9+R19.7+R21+R51+R76.11+R76.8+R80.9+R82.71+Z22.330+Z23+Z29.8+Z33.1+Z34.00+Z34.02+Z34.03+Z34.80+Z34.83+Z34.90+Z34.93+Z36+Z36.3+Z41.8+Z67.91+Z78.9+Z79.4+Z82.49+Z83.3+Z86.19+Z87.440+Z87.891+Z88.0+Z96.41+Z98.890+Z98.891+ACETAMINOPHEN..Acetaminophen.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+AMPICILLIN..Ampicillin.+ASPIRIN..Aspirin.+BECLOMETHASONE.DIPROPIONATE..Beclomethasone.Dipropionate.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+BUPIVACAINE..Bupivacaine.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CEFAZOLIN..Cefazolin.+CEFTRIAXONE..Ceftriaxone.+CEPHALEXIN..Cephalexin.+CETIRIZINE..Cetirizine.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+DOCUSATE.SODIUM..Docusate.Sodium.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+EPINEPHRINE..Epinephrine.+ERYTHROMYCIN..Erythromycin.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.SULFATE..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+FOLIC.ACID..Folic.Acid.+GLUCOSE..Glucose.+Glycol..polyethylene..Polyethylene.Glycols.+HEPARIN..Heparin.+HYDRALAZINE..Hydralazine.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+hydrocortisone.sod.succinate..Hydrocortisone.sodium.succinate.+HYDROMORPHONE..Hydromorphone.+HYDROXYPROGESTERONE.CAPROATE..hydroxyprogesterone.caproate..USP..+INDOMETHACIN..Indomethacin.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LANSOPRAZOLE..lansoprazole.+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORAZEPAM..Lorazepam.+MAGNESIUM.OXIDE..Magnesium.Oxide.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Makena+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+METFORMIN..Metformin.+METOCLOPRAMIDE..Metoclopramide.+MIDAZOLAM..Midazolam.+MORPHINE..Morphine.
                               ,rnd = list(Patient_index=~1),family = family, data = EMR_long_train, lambda=lambda[opt]))
  
    if(class(glm_final)!="try-error"){  
      summary(glm_final)
    
      predictions[[j]] <- predict(glm_final, EMR_long_test, type="response",s=lambda[opt])
      original[[j]]<-EMR_long_test$Term
    }
  }
  save(predictions,original,file=paste0("Results/predictions_all_down_",i,".Rdata"))
}

######################
## Analysis Results ##
######################
library(AUC)
auc_mean<-NULL
for(i in 1:10){
  load(paste0("Results/predictions_all_down_",i,".Rdata"))
  auc<-NULL
  for(j in 1:10){
    print(j)
    if(length(predictions[[j]])!=0){
      auc[j]<-auc(roc(predictions[[j]],factor(original[[j]])))
    }
  }
  auc_mean[i]<-mean(auc,na.rm = T)
}

tiff("Results/AUC_multi_all.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC (labs+diags+meds+demo) - AUC=0.5")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()