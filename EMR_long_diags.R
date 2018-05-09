rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB Diags
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


EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")


#########################
#### Diags data ########
########################
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_diags$Patient_index),"_")
EMR_long_diags$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_diags$Birth_id<- unlist(lapply(splitpop, "[", 2))

##Select all the variables 
EMR_long_diags_data<-EMR_long_diags[,6:(ncol(EMR_long_diags)-2)] ##2321 diags 
rownames(EMR_long_diags_data)<-EMR_long_diags$Sample_ID

table(EMR_long_diags$Term)

matrix<-data.matrix(table(EMR_long_diags$Patient_index,EMR_long_diags$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #155 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #3364 births with Term

###Select the lab test that are complete
num_null<-NULL
for (i in 1:ncol(EMR_long_diags_data)){
  num_null[i]<-dim(table(EMR_long_diags_data[,i]))
}

##None null 

EMR_long_diags_merge<-cbind(EMR_long_diags$Term,EMR_long_diags$WeekOfPregnancy,EMR_long_diags$Patient_index,EMR_long_diags$Individual_id,EMR_long_diags_data)
colnames(EMR_long_diags_merge)[1:4]<-c("Term","WeekOfPregnancy","Pat_birth_id","Patient_id")
EMR_long_diags_merge$Pat_birth_id<-factor(EMR_long_diags_merge$Pat_birth_id)
EMR_long_diags_merge$Term<-factor(EMR_long_diags$Term,levels = c("Term","PTB"))
EMR_long_diags_merge$Patient_id<-factor(EMR_long_diags_merge$Patient_id)

results_2categ<-matrix(NA,ncol(EMR_long_diags_merge),4)
for (i in 5:ncol(EMR_long_diags_merge)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(factor(EMR_long_diags_merge[,i]),ref="0") + WeekOfPregnancy + (1|Pat_birth_id) + (1|Patient_id),
                        data=EMR_long_diags_merge, family=binomial))
  if(class(fm_full)!="try-error"){
    results_2categ[i,1]<-coefficients(summary(fm_full))[2,1] #coef diags
    results_2categ[i,2]<-coefficients(summary(fm_full))[2,4] #p diags
    results_2categ[i,3]<-coefficients(summary(fm_full))[3,1] #coef week
    results_2categ[i,4]<-coefficients(summary(fm_full))[3,4] #p week
    
  }
}
results_2categ<-results_2categ[-c(1:4),]
colnames(results_2categ)<-c("coef_diags","p_diags","coef_week","p_week")
rownames(results_2categ)<-colnames(EMR_long_diags_merge)[-c(1:4)]
write.csv(results_2categ,"results_diags.csv")

##adjust for MT
p_val_long_diags_adj<-p.adjust(results_2categ[,2],method = "fdr")
table(p_val_long_diags_adj<0.05) #129
##Extract significant
id_sign<-match(names(which(p_val_long_diags_adj<0.05)),colnames(EMR_long_diags_merge))
EMR_long_diags_merge[,id_sign]

