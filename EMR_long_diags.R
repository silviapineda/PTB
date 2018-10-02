rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: EMR PTB diags
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


EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")


#########################
#### diags data ########
########################
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_diags$Patient_index),"_")
EMR_long_diags$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_diags$Birth_id<- unlist(lapply(splitpop, "[", 2))
EMR_long_diags$Unique_id<-paste(EMR_long_diags$Patient_index,EMR_long_diags$WeekOfPregnancy,sep="_")
rownames(EMR_long_diags)<-EMR_long_diags$Unique_id


##Select all the variables 
EMR_long_diags_result<-EMR_long_diags[,6:(ncol(EMR_long_diags)-5)] #3157 diags

table(EMR_long_diags$Term) #1928 PTB and 32,291 Term
rownames(EMR_long_diags_result)<-EMR_long_diags$Unique_id

matrix<-data.matrix(table(EMR_long_diags$Patient_index,EMR_long_diags$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #414 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #8648 births with Term

###Filter by nearZeroVar
##Total samples = 34219.
##Uniquecut is the cutoff for the percentage of distinct values out of the number of total samples. 
##FreqCut is the cutoff for the ratio of the most common value to the second most common value. 
n=dim(EMR_long_diags_result)[1]
id_nzv<-nearZeroVar(EMR_long_diags_result,freqCut = (n-10)/10, uniqueCut = 100*(10/n)) 
EMR_long_diags_filter<-EMR_long_diags_result[,-id_nzv] ##855 diags

####Running the univariate longitudinal model
EMR_long_diags$Patient_index<-factor(EMR_long_diags$Patient_index)
EMR_long_diags$Term<-factor(EMR_long_diags$Term,levels = c("Term","PTB"))
EMR_long_diags$Individual_id<-factor(EMR_long_diags$Individual_id)
EMR_long_diags$Outcome<-ifelse(EMR_long_diags$Term=="PTB",1,0)

results_diags<-matrix(NA,ncol(EMR_long_diags_filter),2)
for (i in 1:ncol(EMR_long_diags_filter)){
  print(i)
  fm_full <-  try(glmer(EMR_long_diags$Outcome ~ EMR_long_diags_filter[,i] +
                          (1|EMR_long_diags$Patient_index),
                        family=binomial))
  
  if(class(fm_full)!="try-error"){
    results_diags[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
    results_diags[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
  }
}
colnames(results_diags)<-c("coef_diags","p_diags")
rownames(results_diags)<-colnames(EMR_long_diags_filter)
write.csv(results_diags,"Results/DIAGS/results_diags.csv")

results_diags<-read.csv("Results/DIAGS/results_diags.csv")
results_diags$padj<-p.adjust(results_diags$p_diags)
results_diags$OR_diags<-exp(results_diags$coef_diags)
results_diags_sign<-results_diags[which(results_diags[,3]<0.05),]
id.sign<-match(results_diags_sign$X,colnames(EMR_long_diags_filter))
EMR_long_diags_result_filter_sign<-cbind(EMR_long_diags_filter[,id.sign],EMR_long_diags_filter[,id.sign])
####No significant results

# perc_term<-NULL
# perc_PTB<-NULL
# num_term<-NULL
# num_PTB<-NULL
# for (i in 1:ncol(EMR_long_diags_result_filter_sign)){
#   num_term[i]<-table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,1]
#   num_PTB[i]<-table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,2]
#   perc_term[i]<-(table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,]/ table(EMR_long_diags$Term))[1]
#   perc_PTB[i]<-(table(EMR_long_diags_result_filter_sign[,i],EMR_long_diags$Term)[2,]/ table(EMR_long_diags$Term))[2]
#   test_Term<-EMR_long_diags_result_filter_sign[which(EMR_long_diags$Term=="Term"),i]
#   test_PTB<-EMR_long_diags_result_filter_sign[which(EMR_long_diags$Term=="PTB"),i]
# 
# }
# 
# write.csv(cbind(results_diags_sign,num_term,perc_term,num_PTB,perc_PTB),file="Results/diags/results_diags_sign.csv")
# 
# ####To plot the data
# for(i in 1:ncol(EMR_long_diags_result_filter_sign)){
#   print(i)
#   EMR_long_diags$diags<-EMR_long_diags_result_filter_sign[,i]
#   EMR_long_diags$Term<-factor(EMR_long_diags$Term,levels = c("PTB","Term"))
#   
#  
#   tiff(paste0("Results/diags/",colnames(EMR_long_diags_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
#   p2<-ggplot(EMR_long_diags, aes(x=as.character(WeekOfPregnancy))) +
#     geom_bar(data=EMR_long_diags[EMR_long_diags$Term=="Term",],
#              aes(y=(diags)/length(diags),fill=Term), stat="identity") +
#     geom_bar(data=EMR_long_diags[EMR_long_diags$Term=="PTB",],
#              aes(y=-(diags)/length(diags),fill=Term), stat="identity") +
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
EMR_long_diags_multivariate<-cbind(EMR_long_diags$Outcome,EMR_long_diags$Patient_index,EMR_long_diags_filter)
colnames(EMR_long_diags_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_diags_multivariate,file="Data/EMR_long_diags_multi.Rdata")

