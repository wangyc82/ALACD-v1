# this procedure is for survival analysis for expression level of lncrna, and TCGAbiolinks clinical data
# taken BRCA patients as an example
load("lncrna-count-BRCA-diffLncRNA.RData")

rm(list=setdiff(ls(),c("lncRNA_ldata_BRCA")))
BRCAT2.9_exp<-lncRNA_ldata_BRCA["BRCAT2.9",]

library("TCGAbiolinks")
library("GenomicAlignments")
library("BiocParallel")

clin <- GDCquery_clinic("TCGA-BRCA", type = "clinical")
followUp<-cbind(clin$vital_status,clin$days_to_last_follow_up)
rownames(followUp)<-clin$submitter_id

patient<-substr(colnames(lncRNA_ldata_BRCA),1,12)
followUpA<-followUp[patient,]

event_vec<-rep(1,length(patient))
event_vec[followUpA[,1] %in% "alive"]<-0
Time<-as.numeric(followUpA[,2])


clusMember<-rep(0,length(patient))
clusMember[which(BRCAT2.9_exp>=6)]<-1
clusMember[which(BRCAT2.9_exp<1)]<-2
survData<-data.frame(event_vec,Time,clusMember)
survData<-survData[which(clusMember!=0),]
colnames(survData)<-c("Event","Time","Clus")

library(survival)
pFit <- survfit(Surv(Time,Event) ~ Clus, data = survData);
fit<-survdiff(Surv(Time,Event) ~ Clus, data = survData);
pValue<-1 - pchisq(fit$chisq, length(fit$n) - 1)
plot(pFit,col=c("brown","skyblue"),cex.axis=1.5,xlab="Time (days)", lwd=2, cex=2, ylab="Probability of Overall Survival",cex.lab=1.5);
legend(locator(), paste("P=", round(pValue,digits=4)), col = "black",box.lty=0,cex=1.2);
text(locator(), labels = c("BRCAT2.9 high expression"), col="brown", cex=1.2)
text(locator(), labels = c("BRCAT2.9 low expression"), col="skyblue", cex=1.2)