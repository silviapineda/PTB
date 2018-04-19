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

EMR_long_lab<-read.csv("Data/EMR_LABS_Term_PTB_one_time_point_36_categorical.csv")
EMR_long_diags<-read.csv("Data/EMR_Diagnoses_Term_PTB_longitudinal_36.csv")
EMR_long_meds<-read.csv("Data/EMR_Meds_Term_PTB_longitudinal_36.csv")

####Princiapl component analysis to study data
##Diags
pca <- prcomp(t(EMR_long_diags[,-c(1:5)]))
SPP <- EMR_long_diags[,5]

library("RColorBrewer")
COLOR=brewer.pal(3,"Set1")

pc <- c(1,2)
tiff("pca_plot_EMR_diags.tiff",h=2000,w=2000,res=300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=COLOR[SPP],pch=20,cex=1.7,xlab="PC1",ylab="PC2")
legend("topright", legend=levels(levels.SPP), col=COLOR,pch=20,cex=1.2)
dev.off()

##meds
pca <- prcomp(t(EMR_long_meds[,-c(1:5)]))
SPP <- EMR_long_meds[,5]

library("RColorBrewer")
COLOR=brewer.pal(3,"Set1")

pc <- c(1,2)
tiff("pca_plot_EMR_meds.tiff",h=2000,w=2000,res=300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=COLOR[SPP],pch=20,cex=1.7,xlab="PC1",ylab="PC2")
legend("topright", legend=levels(levels.SPP), col=COLOR,pch=20,cex=1.2)
dev.off()

#########################
#### Diags data ########
########################
EMR_long_diags_data<-EMR_long_diags[,6:ncol(EMR_long_diags)] ##3071 diags
rownames(EMR_long_diags_data)<-EMR_long_diags$Sample_ID
num_diags_PTB<-NULL
num_diags_Term<-NULL
for (i in 1:ncol(EMR_long_diags_data)){
  num_diags_PTB[i]<-table(EMR_long_diags_data[which(EMR_long_diags$Term=="PTB"),i])[2]
  num_diags_Term[i]<-table(EMR_long_diags_data[which(EMR_long_diags$Term=="Term"),i])[2]
}
num_diags_PTB<-ifelse(is.na(num_diags_PTB)==T,0,num_diags_PTB)
num_diags_Term<-ifelse(is.na(num_diags_Term)==T,0,num_diags_Term)

EMR_long_diags_data_PTB<-EMR_long_diags_data[which(EMR_long_diags$Term=="PTB"),which(num_diags_PTB>1)] #567 variables
EMR_long_diags_data_Term<-EMR_long_diags_data[which(EMR_long_diags$Term=="Term"),which(num_diags_Term>1)] #1955 variables

id.common<-match(colnames(EMR_long_diags_data_PTB),colnames(EMR_long_diags_data_Term))
EMR_long_diags_merge<-rbind(EMR_long_diags_data_PTB[,which(is.na(id.common)==F)],EMR_long_diags_data_Term[,na.omit(id.common)]) 
EMR_long_diags_merge<-EMR_long_diags_merge[match(EMR_long_diags$Sample_ID,rownames(EMR_long_diags_merge)),]
##524 variables is the final set with lab test in more than 1 sample per category

p_value_long<-NULL
for (i in 1:ncol(EMR_long_diags_merge)){
  print(i)
  fm_full <-  try(glmer(EMR_long_diags$Term ~ as.numeric(EMR_long_diags_merge[,i]) + EMR_long_diags$WeekOfPregnancy + (1|EMR_long_diags$Patient_ID) ,
                    family=binomial))
  if(class(fm_full)!="try-error"){
    p_value_long[i]<-coefficients(summary(fm_full))[2,4]
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
# num_meds_PTB<-NULL
# num_meds_Term<-NULL
# num_meds<-NULL
# for (i in 1:ncol(EMR_long_meds_data)){
#   num_meds_PTB[i]<-table(EMR_long_meds_data[which(EMR_long_meds$Term=="PTB"),i])[2]
#   num_meds_Term[i]<-table(EMR_long_meds_data[which(EMR_long_meds$Term=="Term"),i])[2]
#   num_meds[i]<-table(EMR_long_meds_data[,i])[2]
# }
# num_meds_PTB<-ifelse(is.na(num_meds_PTB)==T,0,num_meds_PTB)
# num_meds_Term<-ifelse(is.na(num_meds_Term)==T,0,num_meds_Term)
# 
# EMR_long_meds_data_PTB<-EMR_long_meds_data[which(EMR_long_meds$Term=="PTB"),which(num_meds_PTB>1)] #131 variables
# EMR_long_meds_data_Term<-EMR_long_meds_data[which(EMR_long_meds$Term=="Term"),which(num_meds_Term>1)] #338 variables
# 
# id.common<-match(colnames(EMR_long_meds_data_PTB),colnames(EMR_long_meds_data_Term))
# EMR_long_meds_merge<-rbind(EMR_long_meds_data_PTB[,which(is.na(id.common)==F)],EMR_long_meds_data_Term[,na.omit(id.common)]) 
# EMR_long_meds_merge<-EMR_long_meds_merge[match(EMR_long_meds$Sample_ID,rownames(EMR_long_meds_merge)),]
##126 variables is the final set with lab test in more than 1 sample per category

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


