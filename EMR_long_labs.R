rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB labs
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


EMR_long_labs<-read.csv("Data/EMR_LABS_Term_PTB_longitudinal_36_ordered.csv")


################
### LAB TEST ###
###############

##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_labs$Patient_index),"_")
EMR_long_labs$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_labs$Birth_id<- unlist(lapply(splitpop, "[", 2))
EMR_long_labs$Unique_id<-paste(EMR_long_labs$Patient_index,EMR_long_labs$WeekOfPregnancy,sep="_")

##Select all the variables 
EMR_long_labs_order<-EMR_long_labs[,5:417] ##413 lab test 
EMR_long_labs_result<-EMR_long_labs[,418:(ncol(EMR_long_labs)-3)] ##413 lab test
rownames(EMR_long_labs_order)<-EMR_long_labs$Unique_id
rownames(EMR_long_labs_result)<-EMR_long_labs$Unique_id

matrix<-data.matrix(table(EMR_long_labs$Patient_index,EMR_long_labs$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #324 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #5006 births with Term


###Check for the lab test in the order value
num_null_TERM<-NULL
num_null_PTB<-NULL
for (i in 1:ncol(EMR_long_labs_order)){
  num_null_TERM[i]<-table(EMR_long_labs_order[which(EMR_long_labs$Term=="Term"),i])[2]
  num_null_PTB[i]<-table(EMR_long_labs_order[which(EMR_long_labs$Term=="PTB"),i])[2]
}

EMR_long_labs_order_filter_PTB<-EMR_long_labs_order[which(EMR_long_labs$Term=="PTB"),which(num_null_PTB>2)] #129 that has more than 2 samples with result
EMR_long_labs_order_filter_Term<-EMR_long_labs_order[which(EMR_long_labs$Term=="Term"),which(num_null_TERM>2)] #229 that has more than 2 samples with result

id.match<-match(colnames(EMR_long_labs_order_filter_PTB),colnames(EMR_long_labs_order_filter_Term))
PTB_data<-EMR_long_labs_order_filter_PTB[,which(is.na(id.match)==F)]
Term_data<-EMR_long_labs_order_filter_Term[,na.omit(id.match)]
EMR_long_labs_order_filter<-rbind(PTB_data,Term_data)
id.merge<-match(EMR_long_labs$Unique_id,rownames(EMR_long_labs_order_filter))
EMR_long_labs_order_filter<-EMR_long_labs_order_filter[id.merge,]

id<-match(colnames(EMR_long_labs_order_filter),colnames(EMR_long_labs_order))
EMR_long_labs_result_filter<-EMR_long_labs_result[,id]
##125 labs final set

EMR_long_labs_order_filter_num<-EMR_long_labs_order_filter
for(i in 1:ncol(EMR_long_labs_order_filter)){
  EMR_long_labs_order_filter_num[,i]<-ifelse(EMR_long_labs_order_filter[,i]=="Ordered",1,0)
}
  
####Running the longitudinal model
EMR_long_labs$Patient_index<-factor(EMR_long_labs$Patient_index)
EMR_long_labs$Term<-factor(EMR_long_labs$Term,levels = c("Term","PTB"))
EMR_long_labs$Individual_id<-factor(EMR_long_labs$Individual_id)
EMR_long_labs$Outcome<-ifelse(EMR_long_labs$Term=="PTB",1,0)

results_labs<-matrix(NA,ncol(EMR_long_labs_order_filter_num),4)
for (i in 1:ncol(EMR_long_labs_order_filter_num)){
  print(i)
  fm_full <-  try(glmer(EMR_long_labs$Outcome ~ EMR_long_labs_order_filter[,i] + EMR_long_labs_result_filter[,i] +
                           (1|EMR_long_labs$Patient_index),
                         family=binomial))
  
  if(class(fm_full)!="try-error"){
    if(dim(coefficients(summary(fm_full)))[1]==2){
      results_labs[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
      results_labs[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
    } else {
      results_labs[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
      results_labs[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
      results_labs[i,3]<-coefficients(summary(fm_full))[3,1] #coef result
      results_labs[i,4]<-coefficients(summary(fm_full))[3,4] #p result
    }
  }
}
colnames(results_labs)<-c("coef_ordered","p_ordered","coef_result","p_result")
rownames(results_labs)<-colnames(EMR_long_labs_result_filter)
write.csv(results_labs,"Results/LABS/results_labs.csv")

results_labs<-read.csv("Results/LABS/results_labs.csv")
results_labs$padj_order<-p.adjust(results_labs[,3])
results_labs$padj_result<-p.adjust(results_labs[,5])
results_labs$OR_odered<-exp(results_labs$coef_ordered)
results_labs$OR_result<-exp(results_labs$coef_result)

###Check for the significant ones
results_labs_sign<-results_labs[which(results_labs[,3]<0.05 | results_labs[,5]<0.05),]
id.sign<-match(results_labs_sign$X,colnames(EMR_long_labs_result_filter))
EMR_long_labs_result_filter_sign<-EMR_long_labs_result_filter[,id.sign]
EMR_long_labs_order_filter_num_sign<-EMR_long_labs_order_filter_num[,id.sign]
EMR_long_labs_result_filter_sign<-EMR_long_labs_result_filter[,id.sign]
#write.csv(cbind(results_labs_sign,num_term,perc_term,num_PTB,perc_PTB),"results_labs_sign.csv")

perc_term<-NULL
perc_PTB<-NULL
num_term<-NULL
num_PTB<-NULL
mean_Term<-NULL
mean_PTB<-NULL
for (i in 1:ncol(EMR_long_labs_order_filter_num_sign)){
  num_term[i]<-table(EMR_long_labs_order_filter_num_sign[,i],EMR_long_labs$Term)[2,1]
  num_PTB[i]<-table(EMR_long_labs_order_filter_num_sign[,i],EMR_long_labs$Term)[2,2]
  perc_term[i]<-(table(EMR_long_labs_order_filter_num_sign[,i],EMR_long_labs$Term)[2,]/ table(EMR_long_labs$Term))[1]
  perc_PTB[i]<-(table(EMR_long_labs_order_filter_num_sign[,i],EMR_long_labs$Term)[2,]/ table(EMR_long_labs$Term))[2]
  test_Term<-EMR_long_labs_result_filter_sign[which(EMR_long_labs$Term=="Term"),i]
  mean_Term[i]<-mean(test_Term[which(test_Term!=0)])
  test_PTB<-EMR_long_labs_result_filter_sign[which(EMR_long_labs$Term=="PTB"),i]
  mean_PTB[i]<-mean(test_PTB[which(test_PTB!=0)])
}
write.csv(cbind(results_labs_sign,num_term,perc_term,mean_Term,num_PTB,perc_PTB,mean_PTB),file="Results/LABS/results_labs_sign.csv")

####To plot the data
for(i in 1:ncol(EMR_long_labs_result_filter_sign)){
  print(i)
  EMR_long_labs$result<-EMR_long_labs_result_filter_sign[,i]
  EMR_long_labs$ordered<-EMR_long_labs_order_filter_num_sign[,i]
  colnames(EMR_long_labs_result_filter_sign)[i]
  EMR_long_labs$Term<-factor(EMR_long_labs$Term,levels = c("PTB","Term"))

  tiff(paste0("Results/LABS/",colnames(EMR_long_labs_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
    p1<-ggplot(EMR_long_labs[which(EMR_long_labs$result!=0),],aes(x=as.character(WeekOfPregnancy),y=result,fill=Term)) + 
    geom_boxplot()
    print(p1)
  dev.off()

  tiff(paste0("Results/LABS/",colnames(EMR_long_labs_result_filter_sign)[i],"order.tiff"),res=300,w=2000,h=2500)
    p2<-ggplot(EMR_long_labs[which(EMR_long_labs$result>0),], aes(x=as.character(WeekOfPregnancy))) +
    geom_bar(data=EMR_long_labs[EMR_long_labs$Term=="Term",],
           aes(y=(ordered)/length(ordered),fill=Term), stat="identity") +
    geom_bar(data=EMR_long_labs[EMR_long_labs$Term=="PTB",],
           aes(y=-(ordered)/length(ordered),fill=Term), stat="identity") +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip(ylim=c(-0.25,0.25)) +
    scale_y_continuous(breaks=seq(-0.25,0.25,0.125), labels=c(0.25,0.125,0,0.125,0.25)) +
    labs(y="Percentage of ordered", x="Week of pregnancy") +
    ggtitle("                         PTB                                      Term")
    print(p2)
  dev.off()
}
#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
EMR_long_labs_multivariate<-cbind(EMR_long_labs$Outcome,EMR_long_labs$Patient_index,EMR_long_labs_order_filter_num,
                                  EMR_long_labs_result_filter)
colnames(EMR_long_labs_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_labs_multivariate,file="EMR_long_labs_multi_all.Rdata")

##To extract the names of the variables
paste(colnames(EMR_long_labs_multivariate)[1:30],collapse ="+")
paste(colnames(EMR_long_labs_multivariate)[31:ncol(EMR_long_labs_multivariate)],collapse ="+")


#########################
##  Run in the server  ##

load("EMR_long_labs_multi.Rdata")
library(glmmLasso)

##Number of unique patient_index to obtain the train and test set
set.seed(33)
EMR_long_labs_patient<-EMR_long_labs_multivariate[unique(EMR_long_labs_multivariate$Patient_index),c("Term","Patient_index")]
table(EMR_long_labs_patient[,"Term"])
##20% of 318 is 64
##20% of 5012 is 1002
predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:318),63)
  EMR_test_PTB<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:5012),1002)
  EMR_test_Term<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term==0)[id_test_data_Term],]

  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_labs_multivariate$Patient_index)))
  EMR_long_labs_test<-EMR_long_labs_multivariate[id_test_data,] ##3925
  EMR_long_labs_train<-EMR_long_labs_multivariate[-(id_test_data),] ##10491


  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

  BIC_vec<-rep(Inf,length(lambda))
  ## first fit good starting model
  library(MASS);library(nlme)
  #PQL<-glmmPQL(Term~1,random = ~1|Patient_index,family=family,data=EMR_long_labs_train)

  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ Alanine.transaminase_order+Albumin..Serum...Plasma_order+Anion.Gap_order+Atrial.Rate_order+Bile.Acids..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Complement.C3..serum_order+Creat.per.Day..UR_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fibrinogen..Functional_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+
                            Prot.Concentration.UR_order+QTcb_order+RDW_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Ventricular.Rate_order+Alanine.transaminase+Albumin..Serum...Plasma+Anion.Gap+Atrial.Rate+Bile.Acids..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Complement.C3..serum+Creat.per.Day..UR+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Ferritin+Fibrinogen..Functional+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin.A1c+Hemoglobin.A2+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Phosphorus..Serum...Plasma+Platelet.Count+Prot.Concentration.UR+QTcb+RDW+Sedimentation.Rate+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Ventricular.Rate,
                          rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }

  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ Alanine.transaminase_order+Albumin..Serum...Plasma_order+Anion.Gap_order+Atrial.Rate_order+Bile.Acids..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Complement.C3..serum_order+Creat.per.Day..UR_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fibrinogen..Functional_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+
                           Prot.Concentration.UR_order+QTcb_order+RDW_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Ventricular.Rate_order+Alanine.transaminase+Albumin..Serum...Plasma+Anion.Gap+Atrial.Rate+Bile.Acids..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Complement.C3..serum+Creat.per.Day..UR+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Ferritin+Fibrinogen..Functional+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin.A1c+Hemoglobin.A2+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Phosphorus..Serum...Plasma+Platelet.Count+Prot.Concentration.UR+QTcb+RDW+Sedimentation.Rate+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Ventricular.Rate,
                         rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_train, lambda=lambda[opt])
  summary(glm_final)

  predictions[[i]] <- predict(glm_final, EMR_long_labs_test, type="response",s=lambda[opt])
  original[[i]]<-EMR_long_labs_test$Term
}
save(predictions,original,file="predictions_labs.Rdata")


#########################################
### Run LASSO with the whole data set ### 
########################################
set.seed(33)
lambda <- seq(100,0,by=-5)
family = binomial(link = logit)


################## First Simple Method ############################################
## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

BIC_vec<-rep(Inf,length(lambda))
for(j in 1:length(lambda)){
  print(paste("Iteration ", j,sep=""))
  glm1 <- try(glmmLasso(Term ~ Alanine.transaminase_order+Albumin..Serum...Plasma_order+Anion.Gap_order+Atrial.Rate_order+Bile.Acids..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Complement.C3..serum_order+Creat.per.Day..UR_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fibrinogen..Functional_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+
                          Prot.Concentration.UR_order+QTcb_order+RDW_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Ventricular.Rate_order+Alanine.transaminase+Albumin..Serum...Plasma+Anion.Gap+Atrial.Rate+Bile.Acids..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Complement.C3..serum+Creat.per.Day..UR+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Ferritin+Fibrinogen..Functional+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin.A1c+Hemoglobin.A2+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Phosphorus..Serum...Plasma+Platelet.Count+Prot.Concentration.UR+QTcb+RDW+Sedimentation.Rate+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Ventricular.Rate,
                        rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_multivariate, lambda=lambda[j]),silent=TRUE)  
  if(class(glm1)!="try-error"){  
    BIC_vec[j]<-glm1$bic
  }
}

opt<-which.min(BIC_vec)
glm_final <- glmmLasso(Term ~ Alanine.transaminase_order+Albumin..Serum...Plasma_order+Anion.Gap_order+Atrial.Rate_order+Bile.Acids..Total_order+Biophysical.Profile.Score..of.10._order+C.Reactive.Protein_order+Complement.C3..serum_order+Creat.per.Day..UR_order+eGFR.if.African.Amer_order+eGFR.if.non.African.American_order+Eosinophil.Abs.Ct_order+Ferritin_order+Fibrinogen..Functional_order+Gamma.Glutamyl.Transpeptidase_order+Glu.Tol.Post.Glucola_order+Glucose..meter.download_order+Glucose..non.fasting_order+Hematocrit_order+Hemoglobin.A1c_order+Hemoglobin.A2_order+Lipase_order+Lymphocyte.Abs.Cnt_order+Magnesium..Serum...Plasma_order+Monocyte.Abs.Count_order+Neutrophil.Absolute.Count_order+Phosphorus..Serum...Plasma_order+Platelet.Count_order+
                         Prot.Concentration.UR_order+QTcb_order+RDW_order+Sedimentation.Rate_order+Sodium..Serum...Plasma_order+Specific.Gravity_order+Thyroid.Stimulating.Hormone_order+Total.Volume.Collected_order+Ventricular.Rate_order+Alanine.transaminase+Albumin..Serum...Plasma+Anion.Gap+Atrial.Rate+Bile.Acids..Total+Biophysical.Profile.Score..of.10.+C.Reactive.Protein+Complement.C3..serum+Creat.per.Day..UR+eGFR.if.African.Amer+eGFR.if.non.African.American+Eosinophil.Abs.Ct+Ferritin+Fibrinogen..Functional+Gamma.Glutamyl.Transpeptidase+Glu.Tol.Post.Glucola+Glucose..meter.download+Glucose..non.fasting+Hematocrit+Hemoglobin.A1c+Hemoglobin.A2+Lipase+Lymphocyte.Abs.Cnt+Magnesium..Serum...Plasma+Monocyte.Abs.Count+Neutrophil.Absolute.Count+Phosphorus..Serum...Plasma+Platelet.Count+Prot.Concentration.UR+QTcb+RDW+Sedimentation.Rate+Sodium..Serum...Plasma+Specific.Gravity+Thyroid.Stimulating.Hormone+Total.Volume.Collected+Ventricular.Rate,
                       rnd = list(Patient_index=~1),family = family, data = EMR_long_labs_multivariate, lambda=lambda[opt])

save(glm_final,"LASSO_LABS_total.Rdata")


##After running
load("Results/predictions_labs.Rdata")
library(AUC)
auc<-NULL

for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}

conf_matrix<-table(round(predictions[[i]]),factor(original[[i]]))
conf_matrix
library(caret)
sensitivity(conf_matrix)
specificity(conf_matrix)

tiff("Results/AUC_labs.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC labs")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()

