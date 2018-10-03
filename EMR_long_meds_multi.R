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
load("Data/EMR_long_meds_multi.Rdata")

demographics<-read.csv("Data/EMR_patients_patient_race_filtered.csv") 
demographics$Patient_Marital_Status<-gsub("RDP-Dissolved","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("RDP-LG SEP","Unknown/Declined",demographics$Patient_Marital_Status)
demographics$Patient_Marital_Status<-gsub("Widowed","Unknown/Declined",demographics$Patient_Marital_Status)

demographics$Patient_Smoking_Status<-gsub("Heavy Tobacco Smoker","Current Every Day Smoker",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Never Assessed","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Unknown If Ever Smoked","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)
demographics$Patient_Smoking_Status<-gsub("Smoker, Current Status Unknown","Unknown/NeverAssesed",demographics$Patient_Smoking_Status)

EMR_long_meds_multivariate<-EMR_long_meds_multivariate[,-2]
EMR_long_meds_multivariate$ID<-rownames(EMR_long_meds_multivariate)
EMR_meds_demo<-merge(EMR_long_meds_multivariate,demographics,by="Patient_index")
EMR_meds_demo<-subset(EMR_meds_demo,select = -c(ID,X,Term.y))
colnames(EMR_meds_demo)[2]<-"Term"

# ##To extract the names of the variables for the glmmLasso
#paste(colnames(EMR_meds_demo),collapse ="+")


####################################
##  Run glmmLasso for the multi  ##
###################################

##Number of unique patient_index to obtain the train and test set
set.seed(54)
EMR_long_meds_patient<-EMR_meds_demo[match(unique(EMR_meds_demo$Patient_index),EMR_meds_demo$Patient_index)
                                     ,c("Term","Patient_index")]
table(EMR_long_meds_patient[,"Term"])
dim_term<-table(EMR_long_meds_patient[,"Term"])[1] ##2894 unique Term
dim_ptb<-table(EMR_long_meds_patient[,"Term"])[2]  ##218 unique PTB

predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:dim_ptb),0.2*dim_ptb) ##20% of PTB
  EMR_test_PTB<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:dim_term),0.2*dim_term) ##20% of term
  EMR_test_Term<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==0)[id_test_data_Term],]
  
  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_meds_demo$Patient_index)))
  EMR_long_meds_test<-EMR_meds_demo[id_test_data,] 
  EMR_long_meds_train<-EMR_meds_demo[-(id_test_data),] 
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ACETAMINOPHEN..Acetaminophen.+ACETAMINOPHEN.CODEINE..Acetaminophen...Codeine.+ACYCLOVIR..Acyclovir.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+AMPICILLIN..Ampicillin.+ASCORBIC.ACID..Ascorbic.Acid.+ASPIRIN..Aspirin.+AZATHIOPRINE..Azathioprine.+AZITHROMYCIN..Azithromycin.+BECLOMETHASONE.DIPROPIONATE..Beclomethasone.Dipropionate.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+Blood.Sugar..Blood.Glucose.+BUDESONIDE..Budesonide.+BUPIVACAINE..Bupivacaine.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CARBOPROST.TROMETHAMINE..carboprost.tromethamine.+CEFAZOLIN..Cefazolin.+CEFTRIAXONE..Ceftriaxone.+CEPHALEXIN..Cephalexin.+CETIRIZINE..Cetirizine.+CHOLECALCIFEROL..Cholecalciferol.+CITALOPRAM..Citalopram.+CLINDAMYCIN..Clindamycin.+CLOBETASOL..Clobetasol.+CLONAZEPAM..Clonazepam.+CLOTRIMAZOLE..Clotrimazole.+CODEINE.GUAIFENESIN..Codeine...Guaifenesin.+CYCLOBENZAPRINE..cyclobenzaprine.+DEXTROMETHORPHAN.GUAIFENESIN..Dextromethorphan...Guaifenesin.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+Diphenhydramine...Lidocaine+DOCUSATE.SODIUM..Docusate.Sodium.+DOXYLAMINE.SUCCINATE..doxylamine.succinate.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+ENOXAPARIN..Enoxaparin.+EPINEPHRINE..Epinephrine.+ERGOCALCIFEROL..Ergocalciferol.+ERYTHROMYCIN..Erythromycin.+ESCITALOPRAM..Escitalopram.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.GLUCONATE..ferrous.gluconate.+FERROUS.SULFATE..ferrous.sulfate.+ferrous.sulfate.iron..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+FOLIC.ACID..Folic.Acid.+GLUCOSE..Glucose.+GLYBURIDE..Glyburide.+Glycol..polyethylene..Polyethylene.Glycols.+GUAIFENESIN..Guaifenesin.+HEPARIN..Heparin.+HYDRALAZINE..Hydralazine.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+HYDROCORTISONE..Hydrocortisone.+HYDROCORTISONE.ACETATE..hydrocortisone.acetate.+hydrocortisone.sod.succinate..Hydrocortisone.sodium.succinate.+HYDROMORPHONE..Hydromorphone.+HYDROXYCHLOROQUINE..Hydroxychloroquine.+HYDROXYPROGESTERONE.CAPROATE..hydroxyprogesterone.caproate..USP..+hydroxyzine.hcl..Hydroxyzine.Hydrochloride.+INDOMETHACIN..Indomethacin.+INFLIXIMAB..infliximab.+IPRATROPIUM.BROMIDE..Ipratropium.Bromide.+iron..Ferrum.metallicum..Homeopathic.preparation.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LAMOTRIGINE..lamotrigine.+LANSOPRAZOLE..lansoprazole.+LEVOMEFOLATE+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORATADINE..Loratadine.+LORAZEPAM..Lorazepam.+MAGNESIUM..Dietary.Magnesium.+MAGNESIUM.HYDROXIDE..Magnesium.Hydroxide.+MAGNESIUM.OXIDE..Magnesium.Oxide.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Makena+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+MENTHOL..Menthol.+MESALAMINE..mesalamine.+METFORMIN..Metformin.+METHADONE..Methadone.+METHIMAZOLE..Methimazole.+METHYLERGONOVINE..Methylergonovine.+METOCLOPRAMIDE..Metoclopramide.+METOPROLOL.TARTRATE..Metoprolol.Tartrate.+METRONIDAZOLE..Metronidazole.+MICONAZOLE.NITRATE..Miconazole.Nitrate.+MIDAZOLAM..Midazolam.+MISOPROSTOL..Misoprostol.+MONTELUKAST..montelukast.+MORPHINE..Morphine.+Multivitamin..Multivitamin.preparation.+NALOXONE..Naloxone.+NICOTINE..Nicotine.+NIFEDIPINE..Nifedipine.+nifedipine.er+nitrofurantoin.macrocrystals..NITROFURANTOIN..MACROCRYSTALS.+
                            NYSTATIN.TRIAMCINOLONE..Nystatin...Triamcinolone.+OMEPRAZOLE..Omeprazole.+ONDANSETRON..Ondansetron.+ondansetron.hcl..Ondansetron.Hydrochloride.+OSELTAMIVIR..Oseltamivir.+OXYCODONE..Oxycodone.+Oxycodone.Acetaminophen..Acetaminophen...Oxycodone.+PANTOPRAZOLE..pantoprazole.+PHENYLEPHRINE..Phenylephrine.+POTASSIUM.CHLORIDE..Potassium.Chloride.+PREDNISONE..Prednisone.+PROGESTERONE..Progesterone.+Progesterone.Vaginal.Insert+PROMETHAZINE..Promethazine.+Proventil+PYRIDOXINE..Pyridoxine.Drug.Class.+RANITIDINE..Ranitidine.+RHO.D.IMMUNE.GLOBULIN..Rho.D..Immune.Globulin.+Sennosides+SERTRALINE..Sertraline.+SIMETHICONE..Simethicone.+SODIUM.CHLORIDE..Sodium.Chloride.+sodium.citrate.citric.acid..Citric.Acid...sodium.citrate.+SUCRALFATE..Sucralfate.+SULFAMETHOXAZOLE.TRIMETHOPRIM..Trimethoprim.Sulfamethoxazole.Combination.+TERBUTALINE..Terbutaline.+TERCONAZOLE..terconazole.+Tetanus..tetanus.toxoid.vaccine..inactivated.+TRIAMCINOLONE.ACETONIDE..Triamcinolone.Acetonide.+Tuberculin.PPD..Purified.Protein.Derivative.of.Tuberculin.+URSODIOL..Ursodiol.+VALACYCLOVIR..valacyclovir.+ZOLPIDEM..zolpidem.+Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                          ,rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ ACETAMINOPHEN..Acetaminophen.+ACETAMINOPHEN.CODEINE..Acetaminophen...Codeine.+ACYCLOVIR..Acyclovir.+ALBUTEROL.SULFATE..Albuterol.Sulfate.+Aluminum.Hydroxide...Simethicone+AMOXICILLIN..Amoxicillin.+AMPICILLIN..Ampicillin.+ASCORBIC.ACID..Ascorbic.Acid.+ASPIRIN..Aspirin.+AZATHIOPRINE..Azathioprine.+AZITHROMYCIN..Azithromycin.+BECLOMETHASONE.DIPROPIONATE..Beclomethasone.Dipropionate.+BETAMETHASONE.ACETATE..Betamethasone.acetate.+BISACODYL..Bisacodyl.+Blood.Sugar..Blood.Glucose.+BUDESONIDE..Budesonide.+BUPIVACAINE..Bupivacaine.+BUTALBITAL..butalbital.+CALCIUM.CARBONATE..Calcium.Carbonate.+CARBOPROST.TROMETHAMINE..carboprost.tromethamine.+CEFAZOLIN..Cefazolin.+CEFTRIAXONE..Ceftriaxone.+CEPHALEXIN..Cephalexin.+CETIRIZINE..Cetirizine.+CHOLECALCIFEROL..Cholecalciferol.+CITALOPRAM..Citalopram.+CLINDAMYCIN..Clindamycin.+CLOBETASOL..Clobetasol.+CLONAZEPAM..Clonazepam.+CLOTRIMAZOLE..Clotrimazole.+CODEINE.GUAIFENESIN..Codeine...Guaifenesin.+CYCLOBENZAPRINE..cyclobenzaprine.+DEXTROMETHORPHAN.GUAIFENESIN..Dextromethorphan...Guaifenesin.+DEXTROSE..Glucose.+DIPHENHYDRAMINE..Diphenhydramine.+Diphenhydramine...Lidocaine+DOCUSATE.SODIUM..Docusate.Sodium.+DOXYLAMINE.SUCCINATE..doxylamine.succinate.+Edisylate..Prochlorperazine..Prochlorperazine.Edisylate.Salt.+ENOXAPARIN..Enoxaparin.+EPINEPHRINE..Epinephrine.+ERGOCALCIFEROL..Ergocalciferol.+ERYTHROMYCIN..Erythromycin.+ESCITALOPRAM..Escitalopram.+FAMOTIDINE..Famotidine.+FENTANYL..Fentanyl.+FERROUS.GLUCONATE..ferrous.gluconate.+FERROUS.SULFATE..ferrous.sulfate.+ferrous.sulfate.iron..ferrous.sulfate.+FLUCONAZOLE..Fluconazole.+FLUTICASONE..fluticasone.+FOLIC.ACID..Folic.Acid.+GLUCOSE..Glucose.+GLYBURIDE..Glyburide.+Glycol..polyethylene..Polyethylene.Glycols.+GUAIFENESIN..Guaifenesin.+HEPARIN..Heparin.+HYDRALAZINE..Hydralazine.+Hydrocodone.Acetaminophen..Acetaminophen...Hydrocodone.+HYDROCORTISONE..Hydrocortisone.+HYDROCORTISONE.ACETATE..hydrocortisone.acetate.+hydrocortisone.sod.succinate..Hydrocortisone.sodium.succinate.+HYDROMORPHONE..Hydromorphone.+HYDROXYCHLOROQUINE..Hydroxychloroquine.+HYDROXYPROGESTERONE.CAPROATE..hydroxyprogesterone.caproate..USP..+hydroxyzine.hcl..Hydroxyzine.Hydrochloride.+INDOMETHACIN..Indomethacin.+INFLIXIMAB..infliximab.+IPRATROPIUM.BROMIDE..Ipratropium.Bromide.+iron..Ferrum.metallicum..Homeopathic.preparation.+LABETALOL..Labetalol.+lactated.ringers..Lactated.Ringer.s.Solution.+LAMOTRIGINE..lamotrigine.+LANSOPRAZOLE..lansoprazole.+LEVOMEFOLATE+LEVOTHYROXINE..Synthetic.Levothyroxine.+LIDOCAINE..Lidocaine.+LORATADINE..Loratadine.+LORAZEPAM..Lorazepam.+MAGNESIUM..Dietary.Magnesium.+MAGNESIUM.HYDROXIDE..Magnesium.Hydroxide.+MAGNESIUM.OXIDE..Magnesium.Oxide.+MAGNESIUM.SULFATE..Magnesium.Sulfate.+Makena+Maleate..Prochlorperazine..Prochlorperazine.Maleate.+MENTHOL..Menthol.+MESALAMINE..mesalamine.+METFORMIN..Metformin.+METHADONE..Methadone.+METHIMAZOLE..Methimazole.+METHYLERGONOVINE..Methylergonovine.+METOCLOPRAMIDE..Metoclopramide.+METOPROLOL.TARTRATE..Metoprolol.Tartrate.+METRONIDAZOLE..Metronidazole.+MICONAZOLE.NITRATE..Miconazole.Nitrate.+MIDAZOLAM..Midazolam.+MISOPROSTOL..Misoprostol.+MONTELUKAST..montelukast.+MORPHINE..Morphine.+Multivitamin..Multivitamin.preparation.+NALOXONE..Naloxone.+NICOTINE..Nicotine.+NIFEDIPINE..Nifedipine.+nifedipine.er+nitrofurantoin.macrocrystals..NITROFURANTOIN..MACROCRYSTALS.+
                           NYSTATIN.TRIAMCINOLONE..Nystatin...Triamcinolone.+OMEPRAZOLE..Omeprazole.+ONDANSETRON..Ondansetron.+ondansetron.hcl..Ondansetron.Hydrochloride.+OSELTAMIVIR..Oseltamivir.+OXYCODONE..Oxycodone.+Oxycodone.Acetaminophen..Acetaminophen...Oxycodone.+PANTOPRAZOLE..pantoprazole.+PHENYLEPHRINE..Phenylephrine.+POTASSIUM.CHLORIDE..Potassium.Chloride.+PREDNISONE..Prednisone.+PROGESTERONE..Progesterone.+Progesterone.Vaginal.Insert+PROMETHAZINE..Promethazine.+Proventil+PYRIDOXINE..Pyridoxine.Drug.Class.+RANITIDINE..Ranitidine.+RHO.D.IMMUNE.GLOBULIN..Rho.D..Immune.Globulin.+Sennosides+SERTRALINE..Sertraline.+SIMETHICONE..Simethicone.+SODIUM.CHLORIDE..Sodium.Chloride.+sodium.citrate.citric.acid..Citric.Acid...sodium.citrate.+SUCRALFATE..Sucralfate.+SULFAMETHOXAZOLE.TRIMETHOPRIM..Trimethoprim.Sulfamethoxazole.Combination.+TERBUTALINE..Terbutaline.+TERCONAZOLE..terconazole.+Tetanus..tetanus.toxoid.vaccine..inactivated.+TRIAMCINOLONE.ACETONIDE..Triamcinolone.Acetonide.+Tuberculin.PPD..Purified.Protein.Derivative.of.Tuberculin.+URSODIOL..Ursodiol.+VALACYCLOVIR..valacyclovir.+ZOLPIDEM..zolpidem.+Patient_Ethnicity+Patient_Marital_Status+Patient_Smoking_Status+Patient_Age+Patient_Race
                         ,rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_test, lambda=lambda[opt])
  summary(glm_final)
  
  predictions[[i]] <- predict(glm_final, EMR_long_meds_test, type="response",s=lambda[opt])
  original[[i]]<-EMR_long_meds_test$Term
}
save(predictions,original,file="Results/predictions_multi_meds.Rdata")



####
## Analysis Results
####
load("Results/predictions_multi_meds.Rdata")
library(AUC)
auc<-NULL

for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}

conf_matrix<-table(round(predictions[[i]]),factor(original[[i]]))
conf_matrix
sensitivity(conf_matrix)
specificity(conf_matrix)

tiff("Results/AUC_multi_meds.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC meds")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()
