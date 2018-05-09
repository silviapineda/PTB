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
library("RColorBrewer")
library(ggplot2)

working_directory<-"/Users/Pinedasans/PTB/"
setwd(working_directory)


EMR_long_labs<-read.csv("Data/EMR_LABS_Term_PTB_longitudinal_36_ordered.csv")
EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")
EMR_long_meds<-read.csv("Data/EMR_Meds_Term_PTB_longitudinal_36.csv")


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


#########################
#### MEDs data ########
########################
##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_meds$Patient_index),"_")
EMR_long_meds$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_meds$Birth_id<- unlist(lapply(splitpop, "[", 2))

##Select all the variables 
EMR_long_meds_data<-EMR_long_meds[,6:(ncol(EMR_long_meds)-2)] ##427 meds 

table(EMR_long_meds$Term) #426 PTB and 4193 with TERM

matrix<-data.matrix(table(EMR_long_meds$Patient_index,EMR_long_meds$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #218 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #2894 births with Term

###Select the lab test that are complete
num_null<-NULL
for (i in 1:ncol(EMR_long_meds_data)){
  num_null[i]<-dim(table(EMR_long_meds_data[,i]))
}
###None

EMR_long_meds_merge<-cbind(EMR_long_meds$Term,EMR_long_meds$WeekOfPregnancy,EMR_long_meds$Patient_index,
                           EMR_long_meds$Individual_id,EMR_long_meds_data)
colnames(EMR_long_meds_merge)[1:4]<-c("Term","WeekOfPregnancy","Patient_birth_id","Patient_id")
EMR_long_meds_merge$Patient_birth_id<-factor(EMR_long_meds_merge$Patient_birth_id)
EMR_long_meds_merge$Term<-factor(EMR_long_meds_merge$Term,levels = c("Term","PTB"))
EMR_long_meds_merge$Patient_id<-factor(EMR_long_meds_merge$Patient_id)

results_meds<-matrix(NA,ncol(EMR_long_meds_merge),4)
for (i in 5:ncol(EMR_long_meds_merge)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(factor(EMR_long_meds_merge[,i]),ref="0") + WeekOfPregnancy + (1|Patient_birth_id) + 
                          (1|Patient_id),data=EMR_long_meds_merge, family=binomial))
  if(class(fm_full)!="try-error"){
    results_meds[i,1]<-coefficients(summary(fm_full))[2,1] #coef diags
    results_meds[i,2]<-coefficients(summary(fm_full))[2,4] #p diags
    results_meds[i,3]<-coefficients(summary(fm_full))[3,1] #coef week
    results_meds[i,4]<-coefficients(summary(fm_full))[3,4] #p week
    
  }
}
results_meds<-results_meds[-c(1:4),]
colnames(results_meds)<-c("coef_meds","p_meds","coef_week","p_week")
rownames(results_meds)<-colnames(EMR_long_meds_merge)[-c(1:4)]
write.csv(results_meds,"results_meds.csv")

##adjust for MT
p_val_long_meds_adj<-p.adjust(results_meds[,2],method = "fdr")
table(p_val_long_meds_adj<0.05) #6
##Extract significant
id_sign<-match(names(which(p_val_long_meds_adj<0.05)),colnames(EMR_long_meds_merge))
EMR_long_meds_merge_sign<-EMR_long_meds_merge[,c(1:4,id_sign)]
results_meds[which(p_val_long_meds_adj<0.05),]

tiff("Progesterone.Vaginal.Insert_MEDS.tiff",res=300,w=2000,h=2500)
ggplot(EMR_long_meds_merge_sign, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_meds_merge_sign[EMR_long_meds_merge_sign$Term=="Term",], 
           aes(y=(Progesterone.Vaginal.Insert)/length(Progesterone.Vaginal.Insert),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_meds_merge_sign[EMR_long_meds_merge_sign$Term=="PTB",],
           aes(y=-(Progesterone.Vaginal.Insert)/length(Progesterone.Vaginal.Insert),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  coord_flip(ylim=c(-0.05,0.05)) + 
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025), labels=c(0.05,0.025,0,0.025,0.05)) +
  labs(y="Percentage of Progesterone.Vaginal.Insert", x="Week of pregnancy") + 
  ggtitle("                 PTB (10/416)                                      Term (91/4102)")
dev.off()





##Interaction 
p_value_int_med<-NULL
for (i in 4:ncol(EMR_long_meds_merge)){
  print(i)
  fm_full <-  try(glmer(Term ~ EMR_long_meds_merge[,i]*WeekOfPregnancy + (1|Patient_ID) ,data=EMR_long_meds_merge,
                        family=binomial))
  if(class(fm_full)!="try-error"){
    if(dim(coefficients(summary(fm_full)))[1]>3){
      p_value_int_med[i]<-coefficients(summary(fm_full))[4,4]
    }
  }
}

p_val_int_meds_adj<-p.adjust(p_value_int_med,method = "fdr")
table(p_val_int_meds_adj<0.05) #



id_sign<-match(colnames(EMR_long_meds_merge[,which(p_val_int_meds_adj<0.05)]),colnames(EMR_long_meds_merge))

p <- ggplot(fm_full, aes(x = WeekOfPregnancy, y = EMR_long_meds_merge_sign[,i], colour = Term)) +
  geom_point(size=1.2) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),size=0.8) +
  labs(x = "Weeks of Pregnancy",y = "Meds") + theme_bw() + theme_light()

print(p)


################
### LAB TEST ###
###############



##Extract Individual ID for the Patient_index
splitpop <- strsplit(as.character(EMR_long_labs$Patient_index),"_")
EMR_long_labs$Individual_id <- unlist(lapply(splitpop, "[", 1))
EMR_long_labs$Birth_id<- unlist(lapply(splitpop, "[", 2))

##Select all the variables 
EMR_long_labs_data<-EMR_long_labs[,7:(ncol(EMR_long_labs)-2)] ##196 lab test 
rownames(EMR_long_labs_data)<-EMR_long_labs$Sample_ID


matrix<-data.matrix(table(EMR_long_labs$Patient_index,EMR_long_labs$Term))
matrix[which(matrix[,1]!=0 & matrix[,2]!=0),] ##Individuals classify as PTB and Term 
dim(matrix[which(matrix[,1]!=0),]) #106 births with PTB
dim(matrix[which(matrix[,2]!=0),]) #2136 births with Term


###Select the lab test that are complete
num_null<-NULL
for (i in 1:ncol(EMR_long_labs_data)){
  num_null[i]<-dim(table(EMR_long_labs_data[,i]))
}

EMR_long_labs_full<-EMR_long_labs_data[,which(num_null>1)] ##Only 27 lab test are complete

###Obtain the number of categories per lab test
num_categ<-NULL
for (i in 1:ncol(EMR_long_labs_full)){
  num_categ[i]<-dim(table(EMR_long_labs_full[,i]))
}

EMR_long_labs_full_2categ<-EMR_long_labs_full[,which(num_categ==2)]
EMR_long_labs_full_3categ<-EMR_long_labs_full[,which(num_categ==3)]

##For those that only have the categories taken vs taken ab
EMR_long_labs_merge_2categ<-cbind(EMR_long_labs$Term,EMR_long_labs$WeekOfPregnancy,EMR_long_labs$Patient_index,EMR_long_labs$Individual_id,EMR_long_labs_full_2categ)
colnames(EMR_long_labs_merge_2categ)[1:4]<-c("Term","WeekOfPregnancy","Pat_birth_id","Patient_id")
EMR_long_labs_merge_2categ$Pat_birth_id<-factor(EMR_long_labs_merge_2categ$Pat_birth_id)
EMR_long_labs_merge_2categ$Term<-factor(EMR_long_labs_merge_2categ$Term,levels = c("Term","PTB"))
EMR_long_labs_merge_2categ$Patient_id<-factor(EMR_long_labs_merge_2categ$Patient_id)

results_2categ<-matrix(NA,ncol(EMR_long_labs_merge_2categ),4)
for (i in 5:ncol(EMR_long_labs_merge_2categ)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(EMR_long_labs_merge_2categ[,i],ref="Not_taken") + WeekOfPregnancy + (1|Pat_birth_id) + (1|Patient_id),
                        data=EMR_long_labs_merge_2categ, family=binomial))
  if(class(fm_full)!="try-error"){
    results_2categ[i,1]<-coefficients(summary(fm_full))[2,1] #coef not-taken
    results_2categ[i,2]<-coefficients(summary(fm_full))[2,4] #p not taken
    results_2categ[i,3]<-coefficients(summary(fm_full))[3,1] #coef week
    results_2categ[i,4]<-coefficients(summary(fm_full))[3,4] #p week
    
  }
}
results_2categ<-results_2categ[-c(1:4),]
colnames(results_2categ)<-c("coef takenAb","p takenAb","coef week","p week")
rownames(results_2categ)<-colnames(EMR_long_labs_merge_2categ)[-c(1:4)]
write.csv(results_2categ,"results_2categ_labs.csv")

# EMR_long_labs_merge_2categ$Glucose..non.fasting_num<-ifelse(EMR_long_labs_merge_2categ$Glucose..non.fasting=="Not_taken",0,1)
# 
# tiff("LABS_Glucose..non.fasting.tiff",res=300,w=2000,h=2500)
# ggplot(EMR_long_labs_merge_2categ, aes(x=as.character(WeekOfPregnancy))) +
#   geom_bar(data=EMR_long_labs_merge_2categ[EMR_long_labs_merge_2categ$Term=="Term",], 
#            aes(y=(Glucose..non.fasting_num)/length(Glucose..non.fasting_num),fill=Term), stat="identity") +
#   geom_bar(data=EMR_long_labs_merge_2categ[EMR_long_labs_merge_2categ$Term=="PTB",],
#            aes(y=-(Glucose..non.fasting_num)/length(Glucose..non.fasting_num),fill=Term), stat="identity") +
#   geom_hline(yintercept=0, colour="white", lwd=1) +
#   coord_flip(ylim=c(-0.05,0.05)) + 
#   scale_y_continuous(breaks=seq(-0.05,0.05,0.025), labels=c(0.05,0.025,0,0.025,0.05)) +
#   labs(y="Percentage of Abnormal Glucose..non.fasting", x="Week of pregnancy") + 
#   ggtitle("                         PTB                                      Term")
# dev.off()


##For those that have three categories 
EMR_long_labs_merge_3categ<-cbind(EMR_long_labs$Term,EMR_long_labs$WeekOfPregnancy,EMR_long_labs$Patient_index,EMR_long_labs$Individual_id,EMR_long_labs_full_3categ)
colnames(EMR_long_labs_merge_3categ)[1:4]<-c("Term","WeekOfPregnancy","Pat_birth_id","Patient_id")
EMR_long_labs_merge_3categ$Pat_birth_id<-factor(EMR_long_labs_merge_3categ$Pat_birth_id)
EMR_long_labs_merge_3categ$Term<-factor(EMR_long_labs_merge_3categ$Term,levels = c("Term","PTB"))
EMR_long_labs_merge_3categ$Patient_id<-factor(EMR_long_labs_merge_3categ$Patient_id)

##The reference is taken normal
results_3categ<-matrix(NA,ncol(EMR_long_labs_merge_3categ),6)
for (i in 5:ncol(EMR_long_labs_merge_3categ)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(EMR_long_labs_merge_3categ[,i],ref="Taken_normal") + WeekOfPregnancy + (1|Pat_birth_id) +(1|Patient_id),
                        data=EMR_long_labs_merge_3categ, family=binomial))
  if(class(fm_full)!="try-error"){
    results_3categ[i,1]<-coefficients(summary(fm_full))[2,1] #coef not-taken
    results_3categ[i,2]<-coefficients(summary(fm_full))[2,4] #p not taken
    results_3categ[i,3]<-coefficients(summary(fm_full))[3,1] #coef taken-ab
    results_3categ[i,4]<-coefficients(summary(fm_full))[3,4] #p taken-ab
    results_3categ[i,5]<-coefficients(summary(fm_full))[4,1] #coef week
    results_3categ[i,6]<-coefficients(summary(fm_full))[4,4] #p week
    
  }
}
results_3categ<-results_3categ[-c(1:4),]
colnames(results_3categ)<-c("coef notTaken","p notTaken","coef takenAb","p takenAb","coef week","p week")
rownames(results_3categ)<-colnames(EMR_long_labs_merge_3categ)[-c(1:4)]
write.csv(results_3categ,"results_3categ_labs.csv")

EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count<-ifelse(EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count=="Not_taken",1,
                                                          ifelse(EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count=="Taken_normal",0,NA))

tiff("Neutrophil.Absolute.Count_labs_takenNormvsNotTaken.tiff",res=300,w=2000,h=2500)
ggplot(EMR_long_labs_merge_3categ, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="Term",], 
           aes(y=(Neutrophil.Absolute.Count)/length(Neutrophil.Absolute.Count),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="PTB",],
           aes(y=-(Neutrophil.Absolute.Count)/length(Neutrophil.Absolute.Count),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  coord_flip(ylim=c(-0.2,0.2)) + 
  scale_y_continuous(breaks=seq(-0.2,0.2,0.1), labels=c(0.2,0.1,0,0.1,0.2)) +
  labs(y="Percentage of Not Taken vs. Taken Normal (Neutrophil.Absolute.Count)", x="Week of pregnancy") + 
  ggtitle("                         PTB                                      Term")
dev.off()

EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count<-ifelse(EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count=="Taken_normal",0,
                                                             ifelse(EMR_long_labs_merge_3categ$Neutrophil.Absolute.Count=="Taken_abnormal",1,NA))

tiff("Neutrophil.Absolute.Count_labs_takenNormvsTakenAb.tiff",res=300,w=2000,h=2500)
ggplot(EMR_long_labs_merge_3categ, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="Term",], 
           aes(y=(Neutrophil.Absolute.Count)/length(Neutrophil.Absolute.Count),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="PTB",],
           aes(y=-(Neutrophil.Absolute.Count)/length(Neutrophil.Absolute.Count),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  coord_flip(ylim=c(-0.05,0.05)) + 
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025), labels=c(0.05,0.025,0,0.025,0.05)) +
  labs(y="Percentage of Taken Abnormal vs. Taken Normal (Neutrophil.Absolute.Count)", x="Week of pregnancy") + 
  ggtitle("                         PTB                                      Term")
dev.off()
