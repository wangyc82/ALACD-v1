#this procedure is to do the differntial expression analysis

library(limma)

load("lncrna-count-BRCA-diffLncRNA.RData")
# lncrna_ldata_NAcancer is the control data including all tissue aside by tumor
nldata<-cbind(lncRNA_ldata_BRCA,lncrna_ldata_NAcancer)
a<-rep(c("tumor","normal"),c(ncol(lncRNA_ldata_BRCA)),ncol(lncrna_ldata_NAcancer))
design <- model.matrix(~ 0 + factor(substr(a,1,1)))
rownames(design) <- colnames(nldata)
colnames(design) <- c("normal","tumor")
fit<-lmFit(nldata,design)
contrast.matrix <- makeContrasts(tumor-normal,levels=design) # Construct the contrast matrix corresponding to specified contrasts of a set of parameters
fit2 <- contrasts.fit(fit, contrast.matrix) # Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
fit3 <- eBayes(fit2)
Results <- topTable(fit3,number=nrow(fit3))


