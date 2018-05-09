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


matrix<-data.matrix(table(EMR_long_labs$Patient_index,EMR_long_labs$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #324 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #5006 births with Term


###Check for the lab test in the order value
num_null<-NULL
for (i in 1:ncol(EMR_long_labs_order)){
  num_null[i]<-dim(table(EMR_long_labs_order[,i]))
}

###Check for the lab test in the results
num_null<-NULL
for (i in 1:ncol(EMR_long_labs_result)){
 if(dim(table(EMR_long_labs_result[,i]))==2){
   num_null[i]<-table(EMR_long_labs_result[,i])[2]
 } else {
    num_null[i]<-3
  }
} 
EMR_long_labs_order_filter<-EMR_long_labs_order[,which(num_null!=1)] #298 that has more than 1 result
EMR_long_labs_result_filter<-EMR_long_labs_result[,which(num_null!=1)] #298 that has more than 1 result

EMR_long_labs_order_filter_num<-EMR_long_labs_order_filter
for(i in 1:ncol(EMR_long_labs_order_filter)){
  EMR_long_labs_order_filter_num[,i]<-ifelse(EMR_long_labs_order_filter[,1]=="Ordered",1,0)
}
  
####Running the longitudinal model
EMR_long_labs$Patient_index<-factor(EMR_long_labs$Patient_index)
EMR_long_labs$Term<-factor(EMR_long_labs$Term,levels = c("Term","PTB"))
EMR_long_labs$Individual_id<-factor(EMR_long_labs$Individual_id)

###To avoid the error scale the data
#datsc[pvars] <- lapply(datsc[pvars],scale)

results_labs<-matrix(NA,ncol(EMR_long_labs_order_filter_num),4)
for (i in 1:ncol(EMR_long_labs_order_filter_num)){
  print(i)
  fm_full <-  try(glmer(EMR_long_labs$Term ~ EMR_long_labs_order_filter_num[,i] + EMR_long_labs_result_filter[,i] + 
                           (1|EMR_long_labs$Patient_index) + (1|EMR_long_labs$Individual_id),
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
write.csv(results_labs,"results_labs.csv")

results_labs<-read.csv("results_labs.csv")
results_labs_sign<-results_labs[which(results_labs[,3]<0.05 | results_labs[,5]<0.05),]
id.sign<-match(results_labs_sign$X,colnames(EMR_long_labs_result_filter))
EMR_long_labs_result_filter_sign<-EMR_long_labs_result_filter[,id.sign]
EMR_long_labs_order_filter_sign<-EMR_long_labs_order_filter[,id.sign]

perc_term<-NULL
perc_PTB<-NULL
num_term<-NULL
num_PTB<-NULL
for (i in 1:ncol(EMR_long_labs_order_filter_sign)){
  num_term[i]<-table(EMR_long_labs_order_filter_sign[,i],EMR_long_labs$Term)[2,1]
  num_PTB[i]<-table(EMR_long_labs_order_filter_sign[,i],EMR_long_labs$Term)[2,2]
  perc_term[i]<-100*(table(EMR_long_labs_order_filter_sign[,i],EMR_long_labs$Term)[2,]/ table(EMR_long_labs$Term))[1]
  perc_PTB[i]<-100*(table(EMR_long_labs_order_filter_sign[,i],EMR_long_labs$Term)[2,]/ table(EMR_long_labs$Term))[2]
}

####To plot the data
EMR_long_labs$result<-EMR_long_labs_result_filter_sign[,5]
EMR_long_labs$ordered<-ifelse(EMR_long_labs_order_filter_sign[,5]=="Ordered",1,0)

ggplot(EMR_long_labs,aes(x=as.character(WeekOfPregnancy),y=result,fill=Term)) + 
  geom_boxplot()

ggplot(EMR_long_labs, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_labs[EMR_long_labs$Term=="Term",],
           aes(y=(ordered)/length(ordered),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_labs[EMR_long_labs$Term=="PTB",],
           aes(y=-(ordered)/length(ordered),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  coord_flip(ylim=c(-0.25,0.25)) +
  scale_y_continuous(breaks=seq(-0.25,0.25,0.125), labels=c(0.25,0.125,0,0.125,0.25)) +
  labs(y="Percentage of ordered", x="Week of pregnancy") +
  ggtitle("                         PTB                                      Term")



###Multivariate model 
##First: Put everything needed in the same dataframe
EMR_long_labs_multivariate<-cbind(EMR_long_labs$Term,EMR_long_labs$Patient_index,EMR_long_labs_order_filter_sign,EMR_long_labs_result_filter_sign)
colnames(EMR_long_labs_multivariate)[1:2]<-c("Term","Patient_index")
fm_full <-  glmer(EMR_long_labs_multivariate$Term ~  EMR_long_labs_multivariate$Alanine.transaminase_order + EMR_long_labs_multivariate$Chloride..Serum...Plasma_order
                  + EMR_long_labs_multivariate$Complement.C4..serum_order + EMR_long_labs_multivariate$Folate..RBC_order
                  + EMR_long_labs_multivariate$MCHC_order + EMR_long_labs_multivariate$MCV_order + EMR_long_labs_multivariate$Neutrophil.Absolute.Count_order
                  + EMR_long_labs_multivariate$Parvo.Ab.B19.IgM_order + EMR_long_labs_multivariate$RBC.Count_order
                  + EMR_long_labs_multivariate$Rubella.Antibody_order  + EMR_long_labs_multivariate$Ventricular.Rate_order 
                  + EMR_long_labs_multivariate$Vitamin.D..25.Hydroxy_order + EMR_long_labs_multivariate$Alanine.transaminase
                  + EMR_long_labs_multivariate$Chloride..Serum...Plasma
                  + EMR_long_labs_multivariate$Complement.C4..serum + EMR_long_labs_multivariate$Folate..RBC
                  + EMR_long_labs_multivariate$MCHC + EMR_long_labs_multivariate$MCV + EMR_long_labs_multivariate$Neutrophil.Absolute.Count
                  + EMR_long_labs_multivariate$Parvo.Ab.B19.IgM + EMR_long_labs_multivariate$RBC.Count
                  + EMR_long_labs_multivariate$Rubella.Antibody  + EMR_long_labs_multivariate$Ventricular.Rate 
                  + EMR_long_labs_multivariate$Vitamin.D..25.Hydroxy +(1|EMR_long_labs_multivariate$Patient_index),family=binomial)


##Number of unique patient_index to obtain the train and test set
EMR_long_labs_patient<-EMR_long_labs_multivariate[unique(EMR_long_labs_multivariate$Patient_index),c("Term","Patient_index")]
table(EMR_long_labs_patient[,"Term"])

##20% of 309 is 62
##20% of 5021 is 1004
id_test_data_PTB<-sample(c(1:309),62)
EMR_test_PTB<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term=="PTB")[id_test_data_PTB],]
id_test_data_Term<-sample(c(1:5021),1004)
EMR_test_Term<-EMR_long_labs_patient[which(EMR_long_labs_patient$Term=="Term")[id_test_data_Term],]

EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_labs$Patient_index)))
EMR_long_labs_test<-EMR_long_labs[id_test_data,] 
EMR_long_labs_train<-EMR_long_labs[-(id_test_data),] 
EMR_long_labs_order_filter_sign_train<-EMR_long_labs_order_filter_sign[match(EMR_long_labs_train$Unique_id,EMR_long_labs$Unique_id),]
EMR_long_labs_order_filter_sign_test<-EMR_long_labs_order_filter_sign[match(EMR_long_labs_test$Unique_id,EMR_long_labs$Unique_id),]
EMR_long_labs_result_filter_sign_train<-EMR_long_labs_result_filter_sign[match(EMR_long_labs_train$Unique_id,EMR_long_labs$Unique_id),]
EMR_long_labs_result_filter_sign_test<-EMR_long_labs_result_filter_sign[match(EMR_long_labs_test$Unique_id,EMR_long_labs$Unique_id),]

fm_full <-  glmer(EMR_long_labs_train$Term ~ EMR_long_labs_order_filter_sign_train$Activated.Partial.Thromboplastin.Time_order +
                    EMR_long_labs_order_filter_sign_train$Basophil.Abs.Count_order + EMR_long_labs_order_filter_sign_train$Fluid.Volume_order +
                    EMR_long_labs_order_filter_sign_train$MCHC_order + EMR_long_labs_order_filter_sign_train$Platelet.Count_order +
                    EMR_long_labs_order_filter_sign_train$Potassium..Serum...Plasma_order+ EMR_long_labs_order_filter_sign_train$Protein.Creat.Ratio..random_order+
                    EMR_long_labs_order_filter_sign_train$PT_order + EMR_long_labs_order_filter_sign_train$Ventricular.Rate_order + EMR_long_labs_result_filter_sign_train$Activated.Partial.Thromboplastin.Time + 
                    EMR_long_labs_result_filter_sign_train$Basophil.Abs.Count +EMR_long_labs_result_filter_sign_train$Fluid.Volume +
                    EMR_long_labs_result_filter_sign_train$Fluid.Volume +EMR_long_labs_result_filter_sign_train$MCHC + EMR_long_labs_result_filter_sign_train$Platelet.Count +
                    EMR_long_labs_result_filter_sign_train$Potassium..Serum...Plasma + EMR_long_labs_result_filter_sign_train$Potassium..Serum...Plasma +
                    EMR_long_labs_result_filter_sign_train$Protein.Creat.Ratio..random + EMR_long_labs_result_filter_sign_train$PT + EMR_long_labs_result_filter_sign_train$Ventricular.Rate +
                    (1|EMR_long_labs_train$Patient_index),family=binomial)
