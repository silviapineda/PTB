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

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)


EMR_long_meds<-read.csv("Data/EMR_Meds_Term_PTB_longitudinal_36.csv")


#########################
#### meds data ########
########################
##Extract Individual ID for the Patient_index
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_meds$Patient_index),"_")
EMR_long_meds$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_meds$Birth_id<- unlist(lapply(splitpop, "[", 2))
EMR_long_meds$Unique_id<-paste(EMR_long_meds$Patient_index,EMR_long_meds$WeekOfPregnancy,sep="_")

##Select all the variables 
EMR_long_meds_result<-EMR_long_meds[,4:(ncol(EMR_long_meds)-3)] #427 meds

table(EMR_long_meds$Term) #426 PTB and 4,193 Term
rownames(EMR_long_meds_result)<-EMR_long_meds$Unique_id
  
matrix<-data.matrix(table(EMR_long_meds$Patient_index,EMR_long_meds$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #218 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #2894 births with Term

###Select the lab test that are complete
num_null_TERM<-NULL
num_null_PTB<-NULL
for (i in 1:ncol(EMR_long_meds_result)){
  num_null_TERM[i]<-table(EMR_long_meds_result[which(EMR_long_meds$Term=="Term"),i])[2]
  num_null_PTB[i]<-table(EMR_long_meds_result[which(EMR_long_meds$Term=="PTB"),i])[2]
}

EMR_long_meds_result_filter_PTB<-EMR_long_meds_result[which(EMR_long_meds$Term=="PTB"),which(num_null_PTB>2)] #117 that has more than 2 samples with result
EMR_long_meds_result_filter_Term<-EMR_long_meds_result[which(EMR_long_meds$Term=="Term"),which(num_null_TERM>2)] #227 that has more than 2 samples with result

id.match<-match(colnames(EMR_long_meds_result_filter_PTB),colnames(EMR_long_meds_result_filter_Term))
PTB_data<-EMR_long_meds_result_filter_PTB[,which(is.na(id.match)==F)]
Term_data<-EMR_long_meds_result_filter_Term[,na.omit(id.match)]
EMR_long_meds_result_filter<-rbind(PTB_data,Term_data)
id.merge<-match(EMR_long_meds$Unique_id,rownames(EMR_long_meds_result_filter))
EMR_long_meds_result_filter<-EMR_long_meds_result_filter[id.merge,]
##108 meds final set

#@#####Running only the set from Idit
progesterone<-c("HYDROXYPROGESTERONE.CAPROATE..hydroxyprogesterone.caproate..USP..","Makena","PROGESTERONE..Progesterone.","Progesterone.Vaginal.Insert")
delivery_block<-c("ASPIRIN..Aspirin.","INDOMETHACIN..Indomethacin.","MAGNESIUM.SULFATE..Magnesium.Sulfate.",
"NIFEDIPINE..Nifedipine.","TERBUTALINE..Terbutaline.")
##These are not in the filter ones: "nifedipine.er","nifedipine.er","IBUPROFEN..Ibuprofen.",
id_prog<-match(progesterone,colnames(EMR_long_meds_result_filter))
id_delivery<-match(delivery_block,colnames(EMR_long_meds_result_filter))

####Running the longitudinal model
EMR_long_meds$Patient_index<-factor(EMR_long_meds$Patient_index)
EMR_long_meds$Term<-factor(EMR_long_meds$Term,levels = c("Term","PTB"))
EMR_long_meds$Individual_id<-factor(EMR_long_meds$Individual_id)
EMR_long_meds$Outcome<-ifelse(EMR_long_meds$Term=="PTB",1,0)

EMR_long_meds_List<-EMR_long_meds_result_filter[,c(id_prog,id_delivery)]
###Using only Idit List
results_meds<-matrix(NA,ncol(EMR_long_meds_List),2)
for (i in 1:ncol(EMR_long_meds_List)){
  print(i)
  fm_full <-  try(glmer(EMR_long_meds$Outcome ~ EMR_long_meds_List[,i] +
                          (1|EMR_long_meds$Patient_index),
                        family=binomial))
  
  if(class(fm_full)!="try-error"){
    results_meds[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
    results_meds[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
  }
}

####Using all
results_meds<-matrix(NA,ncol(EMR_long_meds_result_filter),2)
for (i in 1:ncol(EMR_long_meds_result_filter)){
  print(i)
  fm_full <-  try(glmer(EMR_long_meds$Outcome ~ EMR_long_meds_result_filter[,i] +
                          (1|EMR_long_meds$Patient_index),
                        family=binomial))
  
  if(class(fm_full)!="try-error"){
    results_meds[i,1]<-coefficients(summary(fm_full))[2,1] #coef ordered
    results_meds[i,2]<-coefficients(summary(fm_full))[2,4] #p ordered
  }
}
colnames(results_meds)<-c("coef_meds","p_meds")
rownames(results_meds)<-colnames(EMR_long_meds_result_filter)
write.csv(results_meds,"Results/MEDS/results_meds.csv")

results_meds<-read.csv("Results/MEDS/results_meds.csv")
results_meds$padj<-p.adjust(results_meds$p_meds)
results_meds$OR_meds<-exp(results_meds$coef_meds)
results_meds_sign<-results_meds[which(results_meds[,3]<0.05),]
id.sign<-match(results_meds_sign$X,colnames(EMR_long_meds_result_filter))
EMR_long_meds_result_filter_sign<-cbind(EMR_long_meds_result_filter[,id.sign],EMR_long_meds_result_filter[,id.sign])

perc_term<-NULL
perc_PTB<-NULL
num_term<-NULL
num_PTB<-NULL
for (i in 1:ncol(EMR_long_meds_result_filter_sign)){
  num_term[i]<-table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,1]
  num_PTB[i]<-table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,2]
  perc_term[i]<-(table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,]/ table(EMR_long_meds$Term))[1]
  perc_PTB[i]<-(table(EMR_long_meds_result_filter_sign[,i],EMR_long_meds$Term)[2,]/ table(EMR_long_meds$Term))[2]
  test_Term<-EMR_long_meds_result_filter_sign[which(EMR_long_meds$Term=="Term"),i]
  test_PTB<-EMR_long_meds_result_filter_sign[which(EMR_long_meds$Term=="PTB"),i]

}

write.csv(cbind(results_meds_sign,num_term,perc_term,num_PTB,perc_PTB),file="Results/MEDS/results_meds_sign.csv")

####To plot the data
for(i in 1:ncol(EMR_long_meds_result_filter_sign)){
  print(i)
  EMR_long_meds$meds<-EMR_long_meds_result_filter_sign[,i]
  EMR_long_meds$Term<-factor(EMR_long_meds$Term,levels = c("PTB","Term"))
  
 
  tiff(paste0("Results/MEDS/",colnames(EMR_long_meds_result_filter_sign)[i],".tiff"),res=300,w=2000,h=2500)
  p2<-ggplot(EMR_long_meds, aes(x=as.character(WeekOfPregnancy))) +
    geom_bar(data=EMR_long_meds[EMR_long_meds$Term=="Term",],
             aes(y=(meds)/length(meds),fill=Term), stat="identity") +
    geom_bar(data=EMR_long_meds[EMR_long_meds$Term=="PTB",],
             aes(y=-(meds)/length(meds),fill=Term), stat="identity") +
    geom_hline(yintercept=0, colour="white", lwd=1) +
    coord_flip(ylim=c(-0.15,0.15)) +
    scale_y_continuous(breaks=seq(-0.15,0.15,0.075), labels=c(0.15,0.075,0,0.075,0.15)) +
    labs(y="Percentage of samples with medication", x="Week of pregnancy") +
    ggtitle("                         PTB                                      Term")
  print(p2)
  dev.off()
}


#############################
### Multivariate model  #####
############################
##First: Put everything needed in the same dataframe
EMR_long_meds_multivariate<-cbind(EMR_long_meds$Outcome,EMR_long_meds$Patient_index,EMR_long_meds_result_filter_sign)
colnames(EMR_long_meds_multivariate)[1:2]<-c("Term","Patient_index")
save(EMR_long_meds_multivariate,file="EMR_long_meds_multi.Rdata")

##To extract the names of the variables
paste(colnames(EMR_long_meds_multivariate),collapse ="+")
fm_full <-  glmer(Term ~ Term+Patient_index+CYCLOBENZAPRINE..cyclobenzaprine.+IRON.SUCROSE..ferric.oxide..saccharated.+Progesterone.Vaginal.Insert+
                    (1|Patient_index),family=binomial,data = EMR_long_meds_multivariate)

#########################
##  Run in the server  ##

load("Data/EMR_long_meds_multi.Rdata")
library(glmmLasso)

##Number of unique patient_index to obtain the train and test set
set.seed(33)
EMR_long_meds_patient<-EMR_long_meds_multivariate[unique(EMR_long_meds_multivariate$Patient_index),c("Term","Patient_index")]
table(EMR_long_meds_patient[,"Term"])
##20% of 300 is 60
##20% of 2812 is 562
predictions<-list() 
original<-list()
for( i in 1:10){
  print(paste("Sample ", i,sep=""))
  id_test_data_PTB<-sample(c(1:300),60)
  EMR_test_PTB<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==1)[id_test_data_PTB],]
  id_test_data_Term<-sample(c(1:2812),562)
  EMR_test_Term<-EMR_long_meds_patient[which(EMR_long_meds_patient$Term==0)[id_test_data_Term],]
  
  EMR_test_data<-rbind(EMR_test_PTB,EMR_test_Term)
  id_test_data<-unlist(lapply(EMR_test_data$Patient_index, function(x) grep(x,EMR_long_meds_multivariate$Patient_index)))
  EMR_long_meds_test<-EMR_long_meds_multivariate[id_test_data,] 
  EMR_long_meds_train<-EMR_long_meds_multivariate[-(id_test_data),] 
  
  
  lambda <- seq(100,0,by=-5)
  family = binomial(link = logit)
  ################## First Simple Method ############################################
  ## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
  
  BIC_vec<-rep(Inf,length(lambda))
  ## first fit good starting model
  library(MASS);library(nlme)
  #PQL<-glmmPQL(Term~1,random = ~1|Patient_index,family=family,data=EMR_long_meds_train)
  
  for(j in 1:length(lambda)){
    print(paste("Iteration ", j,sep=""))
    glm1 <- try(glmmLasso(Term ~ CYCLOBENZAPRINE..cyclobenzaprine.+IRON.SUCROSE..ferric.oxide..saccharated.+Progesterone.Vaginal.Insert,
                          rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_train, lambda=lambda[j]),silent=TRUE)  
    if(class(glm1)!="try-error"){  
      BIC_vec[j]<-glm1$bic
    }
  }
  
  opt<-which.min(BIC_vec)
  glm_final <- glmmLasso(Term ~ CYCLOBENZAPRINE..cyclobenzaprine.+IRON.SUCROSE..ferric.oxide..saccharated.+Progesterone.Vaginal.Insert,
                         rnd = list(Patient_index=~1),family = family, data = EMR_long_meds_train, lambda=lambda[opt])
  summary(glm_final)
  
  #predicting with the training set
  #pred_train <- predict(glm_final, EMR_long_meds_train, type="class")
  #report mean error rate (fraction of incorrect labels)
  #confusion_matrix_train <- table(EMR_long_meds_train$Term, round(pred_train))
  #confusion_matrix_train
  
  #predicting with the test set
  predictions[[i]] <- predict(glm_final, EMR_long_meds_test, type="response",s=lambda[opt])
  #confusion_matrix_test <- table(EMR_long_meds_test$Term, round(pred_test))
  #confusion_matrix_test
  original[[i]]<-EMR_long_meds_test$Term
}
save(predictions,original,file="predictions_meds.Rdata")

##After running
load("Results/predictions_meds.Rdata")
library(AUC)
auc<-NULL
sensitivity<-NULL
for(i in 1:10){
  auc[i]<-auc(roc(predictions[[i]],factor(original[[i]])))
}
conf_matrix<-table(round(predictions[[i]]),factor(original[[i]]))
conf_matrix
library(caret)
sensitivity(conf_matrix)
specificity(conf_matrix)


tiff("Results/AUC_meds.tiff",res=300,w=2000,h=2000)
roc1<-roc(predictions[[1]],factor(original[[1]]))
plot(roc1, col = 1, lty = 2, main = "ROC meds")

for (i in 2:10){
  roci<-roc(predictions[[i]],factor(original[[i]]))
  plot(roci, col = i, lty = 3, add = TRUE)
}
dev.off()



