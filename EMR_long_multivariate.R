rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB
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
library(RColorBrewer)
library(ggplot2)
library(glmmLasso)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)

EMR_long<-read.csv("Data/meds_labs_diagnoses_longitudinal_for_modeling.csv")
table(EMR_long$Term)

##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long$Patient_index),"_")
EMR_long$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long$Birth_id<- unlist(lapply(splitpop, "[", 2))

##Select all the variables 
EMR_long_data<-EMR_long[,6:(ncol(EMR_long)-2)] ##2924 variables 
rownames(EMR_long_data)<-EMR_long_data$Sample_ID

matrix<-data.matrix(table(EMR_long$Patient_index,EMR_long$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #155 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #3364 births with Term

num_null<-NULL
for (i in 1:ncol(EMR_long_data)){
  num_null[i]<-dim(table(EMR_long_data[,i]))
}

EMR_long_data<-EMR_long_data[,which(num_null>1)] ##1827 variables 

###Data frame to run the model
EMR_long_merge<-cbind(EMR_long$Term,EMR_long$WeekOfPregnancy,EMR_long$Patient_index,EMR_long$Individual_id,EMR_long_data)
colnames(EMR_long_merge)[1:4]<-c("Term","WeekOfPregnancy","Pat_birth_id","Patient_id")
EMR_long_merge$Pat_birth_id<-factor(EMR_long_merge$Pat_birth_id)
EMR_long_merge$Term<-factor(EMR_long_merge$Term,levels = c("Term","PTB"))
EMR_long_merge$Patient_id<-factor(EMR_long_merge$Patient_id)

glm.obj <- glmmLasso(fix = Term ~ Neutrophil.Absolute.Count,
                     rnd = list(Pat_birth_id = ~1),family = binomial(), data = EMR_long_merge,
                     lambda=10)
