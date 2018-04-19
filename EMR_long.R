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
setwd("/Users/Pinedasans/PTB/")
library(lattice)
library(lme4)
library("RColorBrewer")
library(ggplot2)

EMR_long_labs<-read.csv("Data/EMR_LABS_Term_PTB_longitudinal_36_categorical.csv")
EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")
EMR_long_meds<-read.csv("Data/EMR_Meds_Term_PTB_longitudinal_36.csv")


#########################
#### Diags data ########
########################
EMR_long_diags_data<-EMR_long_diags[,6:ncol(EMR_long_diags)] ##3066 diags
rownames(EMR_long_diags_data)<-EMR_long_diags$Sample_ID

EMR_long_diags_merge<-cbind(EMR_long_diags$Term,EMR_long_diags$WeekOfPregnancy,EMR_long_diags$Patient_ID,EMR_long_diags_data)
colnames(EMR_long_diags_merge)[1:3]<-c("Term","WeekOfPregnancy","Patient_ID")
EMR_long_diags_merge$Patient_ID<-factor(EMR_long_diags$Patient_ID)
EMR_long_diags_merge$Term<-factor(EMR_long_diags$Term,levels = c("Term","PTB"))

p_value_long_diags<-NULL
for (i in 1:ncol(EMR_long_diags_merge)){
  print(i)
  fm_full <-  try(glmer(Term ~ EMR_long_diags_merge[,i] + WeekOfPregnancy + (1|Patient_ID) ,data=EMR_long_diags_merge,
                    family=binomial))
  if(class(fm_full)!="try-error"){
    p_value_long_diags[i]<-coefficients(summary(fm_full))[2,4]
  }
}

p_val_long_diags_adj<-p.adjust(p_value_long,method = "fdr")
table(p_val_long_diags_adj<0.05) #7

EMR_long_diags[which(EMR_long_diags_merge[,i]==1),1:5]

COLOR=brewer.pal(3,"Set2")
xyplot(EMR_long_diags_merge[,500] ~ EMR_long_diags$WeekOfPregnancy,groups=EMR_long_diags$Term, fill.color = COLOR,pch=20,col = COLOR, 0.5)

#########################
#### MEDs data ########
########################
EMR_long_meds_data<-EMR_long_meds[,6:ncol(EMR_long_meds)] ##457 meds
rownames(EMR_long_meds_data)<-EMR_long_meds$Sample_ID

EMR_long_meds_merge<-cbind(EMR_long_meds$Term,EMR_long_meds$WeekOfPregnancy,EMR_long_meds$Patient_ID,EMR_long_meds_data)
colnames(EMR_long_meds_merge)[1:3]<-c("Term","WeekOfPregnancy","Patient_ID")
EMR_long_meds_merge$Patient_ID<-factor(EMR_long_meds_merge$Patient_ID)
EMR_long_meds_merge$Term<-factor(EMR_long_meds_merge$Term,levels = c("Term","PTB"))

p_value_long<-NULL
for (i in 4:ncol(EMR_long_meds_merge)){
  print(i)
  fm_full <-  try(glmer(Term ~ EMR_long_meds_merge[,i] + WeekOfPregnancy + (1|Patient_ID) ,data=EMR_long_meds_merge,
                        family=binomial))
  if(class(fm_full)!="try-error"){
    p_value_long[i]<-coefficients(summary(fm_full))[2,4]
  }
}

p_val_long_meds_adj<-p.adjust(p_value_long,method = "fdr")
table(p_val_long_meds_adj<0.05) #5 meds significant

id_sign<-match(colnames(EMR_long_meds_merge[,which(p_val_long_meds_adj<0.05)]),colnames(EMR_long_meds_merge))

#"BETAMETHASONE.ACETATE..Betamethasone.acetate." "DIPHENHYDRAMINE..Diphenhydramine."            
#"DOXYCYCLINE.HYCLATE..doxycycline.hyclate."     "MICONAZOLE.NITRATE..Miconazole.Nitrate."      
#"RANITIDINE..Ranitidine."
EMR_long_meds_merge_sign<-EMR_long_meds_merge[,c("Term","WeekOfPregnancy","Patient_ID",colnames(EMR_long_meds_merge)[id_sign])]

tiff("RANITIDINE_meds.tiff",res=300,w=2000,h=2500)
ggplot(EMR_long_meds_merge_sign, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_meds_merge_sign[EMR_long_meds_merge_sign$Term=="Term",], 
           aes(y=(RANITIDINE..Ranitidine.)/length(RANITIDINE..Ranitidine.),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_meds_merge_sign[EMR_long_meds_merge_sign$Term=="PTB",],
           aes(y=-(RANITIDINE..Ranitidine.)/length(RANITIDINE..Ranitidine.),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
          coord_flip(ylim=c(-0.1,0.1)) + 
          scale_y_continuous(breaks=seq(-0.1,0.1,0.05), labels=c(0.1,0.05,0,0.05,0.1)) +
  labs(y="Percentage of RANITIDINE", x="Week of pregnancy") + 
  ggtitle("                         PTB                                      Term")
dev.off()

#extract the model for the significance ones
coef<-matrix(NA,5,4)
for (i in 4:ncol(EMR_long_meds_merge_sign)){
  print(i)
  fm_full <-  try(glmer(Term ~ EMR_long_meds_merge_sign[,i] + WeekOfPregnancy + (1|Patient_ID) ,data=EMR_long_meds_merge_sign,
                        family=binomial))
  
    coef[(i-3),1]<-coefficients(summary(fm_full))[2,1]
    coef[(i-3),2]<-coefficients(summary(fm_full))[2,4]
    coef[(i-3),3]<-coefficients(summary(fm_full))[3,1]
    coef[(i-3),4]<-coefficients(summary(fm_full))[3,4]
}

rownames(coef)<-colnames(EMR_long_meds_merge)[id_sign]
write.csv(coef,"coef_med.csv")

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
EMR_long_labs_data<-EMR_long_labs[,7:ncol(EMR_long_labs)] ##230 labs
rownames(EMR_long_labs_data)<-EMR_long_labs$Sample_ID

num_null<-NULL
for (i in 1:ncol(EMR_long_labs_data)){
  num_null[i]<-dim(table(EMR_long_labs_data[,i]))
}

EMR_long_labs_full<-EMR_long_labs_data[,which(num_null>1)] ##Only 34 lab test are complete
num_categ<-NULL
for (i in 1:ncol(EMR_long_labs_full)){
  num_categ[i]<-dim(table(EMR_long_labs_full[,i]))
}


EMR_long_labs_full_2categ<-EMR_long_labs_full[,which(num_categ==2)]
EMR_long_labs_full_3categ<-EMR_long_labs_full[,which(num_categ==3)]

##For those that only have the categories taken vs taken ab
EMR_long_labs_merge_2categ<-cbind(EMR_long_labs$Term,EMR_long_labs$WeekOfPregnancy,EMR_long_labs$Patient_ID,EMR_long_labs_full_2categ)
colnames(EMR_long_labs_merge_2categ)[1:3]<-c("Term","WeekOfPregnancy","Patient_ID")
EMR_long_labs_merge_2categ$Patient_ID<-factor(EMR_long_labs_merge_2categ$Patient_ID)
EMR_long_labs_merge_2categ$Term<-factor(EMR_long_labs_merge_2categ$Term,levels = c("Term","PTB"))

results_2categ<-matrix(NA,ncol(EMR_long_labs_merge_2categ),4)
for (i in 4:ncol(EMR_long_labs_merge_2categ)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(EMR_long_labs_merge_2categ[,i],ref="Not_taken") + WeekOfPregnancy + (1|Patient_ID),
                        data=EMR_long_labs_merge_2categ, family=binomial))
  if(class(fm_full)!="try-error"){
    results_2categ[i,1]<-coefficients(summary(fm_full))[2,1] #coef not-taken
    results_2categ[i,2]<-coefficients(summary(fm_full))[2,4] #p not taken
    results_2categ[i,3]<-coefficients(summary(fm_full))[3,1] #coef week
    results_2categ[i,4]<-coefficients(summary(fm_full))[3,4] #p week
    
  }
}
results_2categ<-results_2categ[-c(1:3),]
colnames(results_2categ)<-c("coef takenAb","p takenAb","coef week","p week")
rownames(results_2categ)<-colnames(EMR_long_labs_merge_2categ)[-c(1:3)]
write.csv(results_2categ,"results_2categ_labs.csv")

##For those that have three categories 
EMR_long_labs_merge_3categ<-cbind(EMR_long_labs$Term,EMR_long_labs$WeekOfPregnancy,EMR_long_labs$Patient_ID,EMR_long_labs_full_3categ)
colnames(EMR_long_labs_merge_3categ)[1:3]<-c("Term","WeekOfPregnancy","Patient_ID")
EMR_long_labs_merge_3categ$Patient_ID<-factor(EMR_long_labs_merge_3categ$Patient_ID)
EMR_long_labs_merge_3categ$Term<-factor(EMR_long_labs_merge_3categ$Term,levels = c("Term","PTB"))

##The reference is taken normal
results_3categ<-matrix(NA,ncol(EMR_long_labs_merge_3categ),6)
for (i in 4:ncol(EMR_long_labs_merge_3categ)){
  print(i)
  fm_full <-  try(glmer(Term ~ relevel(EMR_long_labs_merge_3categ[,i],ref="Taken_normal") + WeekOfPregnancy + (1|Patient_ID),
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
results_3categ<-results_3categ[-c(1:3),]
colnames(results_3categ)<-c("coef notTaken","p notTaken","coef takenAb","p takenAb","coef week","p week")
rownames(results_3categ)<-colnames(EMR_long_labs_merge_3categ)[-c(1:3)]
write.csv(results_3categ,"results_3categ_labs.csv")


tiff("Lymphocyte.Abs.Cnt_labs.tiff",res=300,w=2000,h=2500)
ggplot(EMR_long_labs_merge_3categ, aes(x=as.character(WeekOfPregnancy))) +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="Term",], 
           aes(y=(Lymphocyte.Abs.Cnt)/length(Lymphocyte.Abs.Cnt),fill=Term), stat="identity") +
  geom_bar(data=EMR_long_labs_merge_3categ[EMR_long_labs_merge_3categ$Term=="PTB",],
           aes(y=-(Lymphocyte.Abs.Cnt)/length(Lymphocyte.Abs.Cnt),fill=Term), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  coord_flip(ylim=c(-0.1,0.1)) + 
  scale_y_continuous(breaks=seq(-0.1,0.1,0.05), labels=c(0.1,0.05,0,0.05,0.1)) +
  labs(y="Percentage of RANITIDINE", x="Week of pregnancy") + 
  ggtitle("                         PTB                                      Term")
dev.off()
