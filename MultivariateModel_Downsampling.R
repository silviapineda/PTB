rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB Diags
###
### CITATION: 
###
### PROCESS: Adjust multivariable model with MEDS+LABS+DIAGS+Demographics
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
library(caret)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)

load("Data/EMR_multi_all.Rdata")

# ##To extract the names of the variables for the glmmLasso
paste(colnames(EMR_labs_diags_meds_filter)[1:100],collapse ="+")
paste(colnames(EMR_labs_diags_meds_filter)[101:ncol(EMR_labs_diags_meds_filter)],collapse ="+")

########################
##### Downsampling ####
######################
set.seed(54)
for (i in 1:10){
  print(paste("down",i))
  ##Number of unique patient_index to obtain the downsampling
  EMR_labs_diags_meds_filter_patient<-EMR_labs_diags_meds_filter[match(unique(EMR_labs_diags_meds_filter$ID),
                                               EMR_labs_diags_meds_filter$ID),c("Term","ID")]
  
  dim_term<-table(EMR_labs_diags_meds_filter_patient[,"Term"])[1]
  dim_ptb<-table(EMR_labs_diags_meds_filter_patient[,"Term"])[2]
  id_down<-sample(c(1:dim_term),dim_ptb)
  EMR_down_term<-EMR_labs_diags_meds_filter_patient[which(EMR_labs_diags_meds_filter_patient$Term==0),][id_down,]
  EMR_down_ptb<-EMR_labs_diags_meds_filter_patient[which(EMR_labs_diags_meds_filter_patient$Term==1),]
  EMR_down<-rbind(EMR_down_term,EMR_down_ptb)
  id_down_data<-unlist(lapply(EMR_down$ID, function(x) grep(x,EMR_labs_diags_meds_filter$ID)))
  EMR_long_diags_multivariate_down<-EMR_labs_diags_meds_filter[id_down_data,]
  
  ##Number of unique patient_index to obtain the train and test set
  EMR_long_diags_multivariate_down_patient<-EMR_long_diags_multivariate_down[match(unique(EMR_long_diags_multivariate_down$ID),
                                                                 EMR_long_diags_multivariate_down$ID),c("Term","ID")]
  dim_term<-table(EMR_labs_diags_meds_filter_patient[,"Term"])[1]
  dim_ptb<-table(EMR_labs_diags_meds_filter_patient[,"Term"])[2]
  
  predictions<-list()
  original<-list()
  set.seed(54)
  for(j in 1:10){
    print(paste("Sample ", j,sep=""))
    id_test_data_PTB<-sample(c(1:dim_ptb),dim_ptb*0.2)
    EMR_test_PTB<-EMR_long_diags_multivariate_down_patient[which(EMR_long_diags_multivariate_down_patient$Term==1)[id_test_data_PTB],]
    id_test_data_Term<-sample(c(1:dim_term),dim_term*0.2)
    EMR_test_Term<-EMR_long_diags_multivariate_down_patient[which(EMR_long_diags_multivariate_down_patient$Term==0)[id_test_data_Term],]
    
    EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
    id_test_data<-unlist(lapply(EMR_test_data$ID, function(x) grep(x,EMR_long_diags_multivariate_down$ID)))
    EMR_long_test<-EMR_long_diags_multivariate_down[id_test_data,] 
    EMR_long_train<-EMR_long_diags_multivariate_down[-(id_test_data),] 
    
    
    lambda <- seq(100,0,by=-5)
    family = binomial(link = logit)
    ################## First Simple Method ############################################
    ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
    
    BIC_vec<-rep(Inf,length(lambda))
    for(k in 1:length(lambda)){
      print(paste("Iteration ", k,sep=""))
      glm1 <- try(glmmLasso(Term ~ X.Int.l.Normaliz.Ratio._order+Activated.Partial.Thromboplastin.Time_order+AFI..cm._order+Alanine.transaminase_order+Albumin..Serum...Plasma_order+Alkaline.Phosphatase_order+Anion.Gap_order+Aspartate.transaminase_order+Atrial.Rate_order+Basophil.Abs.Count_order+Bile.Acids..Total_order+Bilirubin..Total_order+Biophysical.Profile.Score..of.10._order+Calcium..total..Serum...Plasma_order+Calculated.P.Axis_order+Calculated.R.Axis_order+Calculated.T.Axis_order+Carbon.Dioxide..Total_order+Chloride..Serum...Plasma_order+Creatinine_order+Creatinine..random.urine_order+DEEPEST.VERTICAL.POCKET..CM._order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Fetal.Breathing_order+Fetal.Movement_order+Fetal.Tone_order+FHR.Baseline_order+Fibrinogen..Functional_order+Fluid.Volume_order+Free.T4_order+Glucose.Loading.Screen_order+Glucose..120.Mins_order+Glucose..60.Mins_order+Glucose..Fasting..Pregnant..2.hr.GTT_order+Glucose..Fasting..Pregnant..3.hr.GTT_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin_order+Hemoglobin.A1c_order+Hours.Collected_order+Lymphocyte.Abs.Cnt_order+MCH_order+MCHC_order+MCV_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Nonstress.Test_order+OGT.Glucola.Dose_order+P.R.Interval_order+pH..UA_order+Platelet.Count_order+POCT.pH..UA_order+POCT.Spec.Grav..UA_order+Potassium..Serum...Plasma_order+Prot.Concentration.UR_order+Prot.Total.per.Day.UR_order+Protein.Creat.Ratio..random_order+PT_order+QRS.Duration_order+QT.Interval_order+QTcb_order+RBC.Count_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Urea.Nitrogen..Serum...Plasma_order+Ventricular.Rate_order+WBC.Count_order+X.Int.l.Normaliz.Ratio.+Activated.Partial.Thromboplastin.Time+AFI..cm.+Alanine.transaminase+Alkaline.Phosphatase+Anion.Gap+Aspartate.transaminase+Atrial.Rate+Basophil.Abs.Count+Bilirubin..Total+Biophysical.Profile.Score..of.10.+Calcium..total..Serum...Plasma+Calculated.P.Axis+Calculated.R.Axis+Calculated.T.Axis+Carbon.Dioxide..Total+Chloride..Serum...Plasma+Creatinine+Creatinine..random.urine+DEEPEST.VERTICAL.POCKET..CM.+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Fetal.Breathing+Fetal.Movement+Fetal.Tone+FHR.Baseline+Fibrinogen..Functional+
                              Fluid.Volume+Free.T4+Glucose.Loading.Screen+Glucose..120.Mins+Glucose..60.Mins+Glucose..Fasting..Pregnant..2.hr.GTT+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin+Hemoglobin.A1c+Hours.Collected+Lymphocyte.Abs.Cnt+MCH+MCHC+MCV+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Nonstress.Test+OGT.Glucola.Dose+P.R.Interval+pH..UA+Platelet.Count+POCT.pH..UA+POCT.Spec.Grav..UA+Potassium..Serum...Plasma+Prot.Concentration.UR+Prot.Total.per.Day.UR+Protein.Creat.Ratio..random+PT+QRS.Duration+QT.Interval+QTcb+RBC.Count+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Urea.Nitrogen..Serum...Plasma+Ventricular.Rate+WBC.Count+D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.3+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O47.00+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891+ACETAMINOPHEN..Acetaminophen.+ACYCLOVIR..Acyclovir.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+ASCORBIC.ACID..Ascorbic.Acid.+ASPIRIN..Aspirin.+AZITHROMYCIN..Azithromycin.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+Blood.Sugar..Blood.Glucose.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CLOTRIMAZOLE..Clotrimazole.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+DOCUSATE.SODIUM..Docusate.Sodium.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.SULFATE..ferrous.sulfate.+ferrous.sulfate.iron..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+GLUCOSE..Glucose.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORAZEPAM..Lorazepam.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+METFORMIN..Metformin.+METRONIDAZOLE..Metronidazole.+NIFEDIPINE..Nifedipine.+nitrofurantoin.macrocrystals..NITROFURANTOIN..MACROCRYSTALS.+ondansetron.hcl..Ondansetron.Hydrochloride.+PREDNISONE..Prednisone.+RANITIDINE..Ranitidine.+RHO.D.IMMUNE.GLOBULIN..Rho.D..Immune.Globulin.+Sennosides+SERTRALINE..Sertraline.+SIMETHICONE..Simethicone.+SODIUM.CHLORIDE..Sodium.Chloride.+Tetanus..tetanus.toxoid.vaccine..inactivated.+TRIAMCINOLONE.ACETONIDE..Triamcinolone.Acetonide.+Tuberculin.PPD..Purified.Protein.Derivative.of.Tuberculin.
                            ,rnd = list(ID=~1),family = family, data = EMR_long_train, lambda=lambda[j]),silent=TRUE)  
      if(class(glm1)!="try-error"){  
        BIC_vec[j]<-glm1$bic
      }
    }
    
    opt<-which.min(BIC_vec)
    glm_final <- try(glmmLasso(Term ~ D64.9+E03.9+E07.9+E66.9+F32.9+F41.9+I10+IMO0001+IMO0002+O09.293+O09.299+O09.513+O09.519+O09.523+O09.529+O09.813+O09.819+O09.893+O09.899+O09.90+O09.93+O10.019+O10.913+O13.3+O13.9+O24.410+O24.414+O24.419+O24.919+O26.843+O26.849+O26.893+O28.9+O32.1XX0+O34.219+O35.8XX0+O35.9XX0+O36.0990+O36.5990+O36.8990+O40.9XX0+O43.93+O44.00+O44.03+O47.00+O69.89X0+O99.019+O99.119+O99.210+O99.213+O99.280+O99.283+O99.810+O99.89+R76.11+Z11.1+Z23+Z30.09+Z33.1+Z34.00+Z34.03+Z34.80+Z34.82+Z34.83+Z34.90+Z34.92+Z34.93+Z36+Z36.3+Z36.9+Z39.1+Z79.4+Z98.890+Z98.891
                               ,rnd = list(ID=~1),family = family, data = EMR_long_diags_train, lambda=lambda[opt]))
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
