source("D:/Zhou, Wenru/Shared Code/temp_table1.R")
library(Hmisc)
library(rlist)
library(DESeq2)
library(DEFormats)
library(EDASeq)
library(edgeR)
library(MASS)
library(gnlm)
library(NBPSeq)
library(dplyr)
library(jrnoldmisc)
library(scales)
library(gnlm)

ipf<-readRDS("D:/Course/BIOS7695/Final projects/data/ipf_data_for_wenru.RDS", refhook = NULL)
demographics<-ipf@colData@listData

demo.dat<-data.frame(matrix(unlist(demographics), nrow=340, byrow=F),stringsAsFactors=FALSE,row.names=ipf@colData@rownames)
colnames(demo.dat)<-rownames(summary(ipf@colData@listData))

demo.dat[demo.dat==""]<-NA

demo.dat$female<-as.factor(demo.dat$female)

demo.dat$case<-as.factor(demo.dat$case)
levels(demo.dat$case)<-c("Control","Case")

demo.dat$PC1_MEGA<-as.numeric(demo.dat$PC1_MEGA)
demo.dat$PC2_MEGA<-as.numeric(demo.dat$PC2_MEGA)
demo.dat$PC3_MEGA<-as.numeric(demo.dat$PC3_MEGA)

demo.dat$Ever.Smoked[demo.dat$Ever.Smoked=="Y"]<-1
demo.dat$Ever.Smoked[demo.dat$Ever.Smoked=="N"]<-0
demo.dat$Ever.Smoked<-as.factor(demo.dat$Ever.Smoked)

demo.dat$Current.smoker[demo.dat$Current.smoker=="Y"]<-1
demo.dat$Current.smoker[demo.dat$Current.smoker=="N"]<-0
demo.dat$Current.smoker[demo.dat$Current.smoker==" "]<-NA
demo.dat$Current.smoker<-as.factor(demo.dat$Current.smoker)

demo.dat$race[demo.dat$race=="black"]<-"Black or African American"
demo.dat$race[demo.dat$race=="white"]<-"White"
demo.dat$race[demo.dat$race=="asian"]<-"Asian"
demo.dat$race[demo.dat$race=="investigated-unknown"|demo.dat$race=="multiple-unspecified"|demo.dat$race=="Multiple Unspecified"|demo.dat$race=="unknown"]<-NA

demo.dat$ethnicity[demo.dat$ethnicity=="hispanic"|demo.dat$ethnicity=="hispanic or latino"]<-"Hispanic or Latino"
demo.dat$ethnicity[demo.dat$ethnicity=="non-hispanic"|demo.dat$ethnicity=="non-hispanic or latino"|demo.dat$ethnicity=="not hispanic or latino"|demo.dat$ethnicity=="Not Hispanic or Latino"]<-"Not Hispanic or Latino"

label(demo.dat$female)<-"Female"
label(demo.dat$case)<-"Case"
label(demo.dat$Ever.Smoked)<-"Ever Smoked"
label(demo.dat$Current.smoker)<-"Current smoker"
label(demo.dat$PC1_MEGA)<-"PC1 MEGA"
label(demo.dat$PC2_MEGA)<-"PC2 MEGA"
label(demo.dat$PC3_MEGA)<-"PC3 MEGA"
label(demo.dat$race)<-"Race"
label(demo.dat$ethnicity)<-"Ethnicity"
label(demo.dat$replaceable)<-"Replaceable"


#Test histrogram
hist(demo.dat$PC1_MEGA)
hist(demo.dat$PC2_MEGA)
hist(demo.dat$PC3_MEGA)
colnames(demo.dat)


table(demo.dat$female)
table(demo.dat$case)
table(demo.dat$Ever.Smoked)
table(demo.dat$Current.smoker)
table(demo.dat$race)
table(demo.dat$ethnicity)
table(demo.dat$replaceable)

#tab.dem<-final_table(demo.dat,c("female","Ever.Smoked","Current.smoker","PC1_MEGA","PC2_MEGA","PC3_MEGA","race","ethnicity","replaceable"),demo.dat$case,margin=2,single=F,2,col.names=T, summary.stat='median')
demo.dat0<-demo.dat
tab.dem0<-final_table(demo.dat0,c("female","Ever.Smoked","Current.smoker","PC1_MEGA","PC2_MEGA","PC3_MEGA","race","ethnicity","replaceable"),demo.dat$case,margin=2,single=F,ron=3,col.names=T, summary.stat='both')

missing.dem<-missing_table(demo.dat0,c("female","Ever.Smoked","Current.smoker","PC1_MEGA","PC2_MEGA","PC3_MEGA","race","ethnicity","replaceable"),demo.dat$case)



#Filter the counts

#Exclude gene rows with counts>=10
detect.zero<-rowSums(matrix(counts(ipf) %in% seq(from=0,to=10,by=1),nrow=nrow(counts(ipf))))
table(detect.zero==0)
#KEEP detect.zero==0
filter10<-counts(ipf)[detect.zero==0,]
test<-rowSums(matrix(filter10 %in% seq(from=0,to=10,by=1),nrow=nrow(filter10)))
table(test==0)

#Exclude participants with missing ever smoke
keeplist<-rownames(demo.dat)[is.na(demo.dat$Ever.Smoked)==FALSE]
demo.dat<-demo.dat[is.na(demo.dat$Ever.Smoked)==FALSE,]
filter10<-filter10[,colnames(filter10) %in% keeplist]
counts<-filter10

#Calculate the demo table again
tab.dem<-final_table(demo.dat,c("female","Ever.Smoked","Current.smoker","PC1_MEGA","PC2_MEGA","PC3_MEGA","race","ethnicity","replaceable"),demo.dat$case,margin=2,single=F,ron=3,col.names=T, summary.stat='both')
