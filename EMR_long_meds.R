rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB meds
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
library(caret)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)


EMR_long_meds<-read.csv("Data/EMR_Meds_Term_PTB_longitudinal_36.csv")


#########################
#### meds data ########
########################
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_meds$Patient_index),"_")
EMR_long_meds$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_meds$Birth_id<- unlist(lapply(splitpop, "[", 2))
EMR_long_meds$Unique_id<-paste(EMR_long_meds$Patient_index,EMR_long_meds$WeekOfPregnancy,sep="_")
rownames(EMR_long_meds)<-EMR_long_meds$Unique_id


##Select all the variables 
EMR_long_meds_result<-EMR_long_meds[,4:(ncol(EMR_long_meds)-3)] #427 meds

table(EMR_long_meds_result$Term) #426 PTB and 4,193 Term
rownames(EMR_long_meds_result)<-EMR_long_meds$Unique_id
  
matrix<-data.matrix(table(EMR_long_meds_result$Patient_index,EMR_long_meds_result$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #218 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #2894 births with Term

###Filter by nearZeroVar
id_nzv<-nearZeroVar(EMR_long_meds_result,freqCut = 99/1,uniqueCut = 1)
EMR_long_meds_filter<-EMR_long_meds_result[,-id_nzv] ##50 meds

####Running the univariate longitudinal model
EMR_long_meds$Patient_index<-factor(EMR_long_meds$Patient_index)
EMR_long_meds$Term<-factor(EMR_long_meds$Term,levels = c("Term","PTB"))
EMR_long_meds$Individual_id<-factor(EMR_long_meds$Individual_id)
EMR_long_meds$Outcome<-ifelse(EMR_long_meds$Term=="PTB",1,0)

results_meds<-matrix(NA,ncol(EMR_long_meds_filter),2)
for (i in 1:ncol(EMR_long_meds_filter)){
  print(i)
  fm_full <-  try(glmer(EMR_long_meds$Outcome ~ EMR_long_meds_filter[,i] +
                          (1|EMR_long_meds$Patient_index),
                        family=binomial))
  
  if(class(fm_full)!="try-error"){
    results_meds[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
    results_meds[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
  }
}
colnames(results_meds)<-c("coef_meds","p_meds")
rownames(results_meds)<-colnames(EMR_long_meds_filter)
write.csv(results_meds,"Results/MEDS/results_meds.csv")

results_meds<-read.csv("Results/MEDS/results_meds.csv")
results_meds$padj<-p.adjust(results_meds$p_meds)
results_meds$OR_meds<-exp(results_meds$coef_meds)
results_meds_sign<-results_meds[which(results_meds[,3]<0.05),]
id.sign<-match(results_meds_sign$X,colnames(EMR_long_meds_filter))
EMR_long_meds_result_filter_sign<-cbind(EMR_long_meds_filter[,id.sign],EMR_long_meds_filter[,id.sign])
####No significant results

# perc_term<-NULL
# perc_PTB<-NULL
# num_term<-NULL
# num_PTB<-NULL
# for (i in 1:ncol(EMR_long_meds_result_filter_sign)){
#   num_term[i]<-table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,1]
#   num_PTB[i]<-table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,2]
#   perc_term[i]<-(table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,]/ table(EMR_long_meds$Term))[1]
#   perc_PTB[i]<-(table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,]/ table(EMR_long_meds$Term))[2]
#   test_Term<-EMR_long_meds_result_filter_sign[which(EMR_long_meds$Term=="Term"),i]
#   test_PTB<-EMR_long_meds_result_filter_sign[which(EMR_long_meds$Term=="PTB"),i]
# 
# }
# 
# write.csv(cbind(results_meds_sign,num_term,perc_term,num_PTB,perc_PTB),file="Results/MEDS/results_meds_sign.csv")
# 
# ####To plot the data
# for(i in 1:ncol(EMR_long_meds_result_filter_sign)){
#   print(i)
#   EMR_long_meds$meds<-EMR_long_meds_result_filter_sign[,i]
#   EMR_long_meds$Term<-factor(EMR_long_meds$Term,levels = c("PTB","Term"))
#   
#  
#   tiff(paste0("Results/MEDS/",colnames(EMR_long_meds_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
#   p2<-ggplot(EMR_long_meds, aes(x=as.character(WeekOfPregnancy))) +
#     geom_bar(data=EMR_long_meds[EMR_long_meds$Term=="Term",],
#              aes(y=(meds)/length(meds),fill=Term), stat="identity") +
#     geom_bar(data=EMR_long_meds[EMR_long_meds$Term=="PTB",],
#              aes(y=-(meds)/length(meds),fill=Term), stat="identity") +
#     geom_hline(yintercept=0, colour="white", lwd=1) +
#     coord_flip(ylim=c(-0.15,0.15)) +
#     scale_y_continuous(breaks=seq(-0.15,0.15,0.075), labels=c(0.15,0.075,0,0.075,0.15)) +
#     labs(y="Percentage of samples with medication", x="Week of pregnancy") +
#     ggtitle("                         PTB                                      Term")
#   print(p2)
#   dev.off()
# }


#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
EMR_long_meds_multivariate<-cbind(EMR_long_meds$Outcome,EMR_long_meds$Patient_index,EMR_long_meds_filter)
colnames(EMR_long_meds_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_meds_multivariate,file="Data/EMR_long_meds_multi.Rdata")

