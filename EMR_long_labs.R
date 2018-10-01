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
### DESCRIP: Analysis of EMR data for PTB lab test univariately
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

###Filter by nearZeroVar
##Total samples = 13214 At least 10 values different. freqCut = 13214-19[10+9]=13195, uniqueCut = (10/13214)*100 
n=dim(EMR_long_labs_order)[1]
id_nzv_order<-nearZeroVar(EMR_long_labs_order,freqCut = n-19,uniqueCut = 100*(10/n))
id_nzv_result<-nearZeroVar(EMR_long_labs_result,freqCut = n-19,uniqueCut = 100*(10/n))
match(id_nzv_order,id_nzv_result) ##All the ones in order are included in results, so taken the result
###Build filtered mayrix
EMR_long_labs_order_filter<-EMR_long_labs_order[,-id_nzv_result] 
EMR_long_labs_result_filter<-EMR_long_labs_result[,-id_nzv_result]

##189 labs final set
###Converting to number
EMR_long_labs_order_filter_num<-EMR_long_labs_order_filter
for(i in 1:ncol(EMR_long_labs_order_filter)){
  EMR_long_labs_order_filter_num[,i]<-ifelse(EMR_long_labs_order_filter[,i]=="Ordered",1,0)
}
  
####Running the univariate longitudinal model
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
write.csv(results_labs,"Results/results_labs.csv")

results_labs<-read.csv("Results/results_labs.csv")
results_labs$padj_order<-p.adjust(results_labs[,3])
results_labs$padj_result<-p.adjust(results_labs[,5])
results_labs$OR_odered<-exp(results_labs$coef_ordered)
results_labs$OR_result<-exp(results_labs$coef_result)

###Check for the significant ones
results_labs_sign<-results_labs[which(results_labs[,6]<0.05 | results_labs[,7]<0.05),]
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
write.csv(cbind(results_labs_sign,num_term,perc_term,mean_Term,num_PTB,perc_PTB,mean_PTB),file="Results/results_labs_sign.csv")

####To plot the data
for(i in 1:ncol(EMR_long_labs_result_filter_sign)){
  print(i)
  EMR_long_labs$result<-EMR_long_labs_result_filter_sign[,i]
  EMR_long_labs$ordered<-EMR_long_labs_order_filter_num_sign[,i]
  colnames(EMR_long_labs_result_filter_sign)[i]
  EMR_long_labs$Term<-factor(EMR_long_labs$Term,levels = c("PTB","Term"))

  tiff(paste0("Results/",colnames(EMR_long_labs_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
    p1<-ggplot(EMR_long_labs[which(EMR_long_labs$result!=0),],aes(x=as.character(WeekOfPregnancy),y=result,fill=Term)) + 
    geom_boxplot()
    print(p1)
  dev.off()

  tiff(paste0("Results/",colnames(EMR_long_labs_result_filter_sign)[i],"order.tiff"),res=300,w=2000,h=2500)
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

###################################################
### Prepare data for the Multivariate model  #####
##################################################
##First: Put everything needed in the same dataframe
EMR_long_labs_multivariate<-cbind(EMR_long_labs$Outcome,EMR_long_labs$Patient_index,EMR_long_labs_order_filter_num,
                                  EMR_long_labs_result_filter)
colnames(EMR_long_labs_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_labs_multivariate,file="Data/EMR_long_labs_multi.Rdata")










