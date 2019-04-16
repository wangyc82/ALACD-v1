# this procedure is to generate intital relationships between lncRNA and mRNA based on their correlation analysis
# taken BRCA TCGA patients as an example
load("BRCA-example.RData")

#lncRNA_ldata_BRCA is expression profiles of lncRNAs across BRCA tumors in TCGA
#GExp_match is expression profiles of mRNAs across same BRCA tumors in TCGA

#creat coexpression matrix
M<-rbind(lncRNA_ldata_BRCA,GExp_match)
CC<-cor(M)
corMat<-CC[1:nrow(lncRNA_ldata_BRCA),(nrow(lncRNA_ldata_BRCA)+1):nrow(CC));
corVec<-as.vector(corMat)
corVec<-abs(corVec)

hl<-hist(corVec)
plot(hl$count, log="y", bty="n",type='h', lwd=10, lend=2,cex.axis=1.5,cex.lab=1.5,col=c("azure4",rep("azure2",22),rep("mistyrose",10)),xlab="Correaltion Coefficient",main="Hisogram of correlation coefficient",ylab="Frequency",xaxt="n")
axis(1, 1:34, hl$breaks,cex.axis=1.2)
#hist(corVec[corVec>=0.2],cex.axis=1.5,cex.lab=1.5,col=c(rep("azure2",5),rep("mistyrose",11)),xlab="Correaltion Coefficient",main="Hisogram of correlation coefficient")

# scale the expression data
GExp_match_scale<-t(apply(GExp_match,1,function(x) 2*(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))))
lncRNA_match_scale<-t(apply(lncRNA_exp_BRCA_match,1,function(x) 2*(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))))

A<-which(abs(corMat)>=0.7,arr.ind = TRUE) # return positives
rowlabP<-A[,1]
collabP<-A[,2]

B<-which(abs(corMat)<=0.000002,arr.ind = TRUE) # return negative candidates
rowlabU<-B[,1]
collabU<-B[,2]

Xp<-cbind(lncRNA_match_scale[rowlabP,],GExp_match_scale[collabP,])
Xn<-cbind(lncRNA_match_scale[rowlabU,],GExp_match_scale[collabU,])
X<-rbind(Xp,Xn)
Ylab<-c(rep("a",nrow(A)),rep("b",nrow(B)))

source('crossvalSVM.R')
crossV<-crossvalSVM(X,Ylab,10,10,0.001)
source('~getperf.R')
perf<-getperf(crossV[[1]],crossV[[2]])


#for prediction
library(e1071)
library(pracma)

model<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.001,class.weights=c(a=1, b=1),probability=TRUE)
preP<-predict(model, Xp,probability=TRUE)
predictedYP<-attr(preP,"probabilities")[,1]

C<-which(abs(corMat)>=0.5 & abs(corMat)<0.7,arr.ind = TRUE) # return abs(cor) fall between 0.3 and 0.5
rowlabT<-C[,1]
collabT<-C[,2]
p<-randperm(length(rowlabT),10000)
pre1<-predict(model,cbind(lncRNA_match_scale[rowlabT[p],],GExp_match_scale[collabT[p],]),probability=TRUE)
predictedY1<-attr(pre1,"probabilities")[,1]

C1<-which(abs(corMat)>=0.1 & abs(corMat)<0.2,arr.ind = TRUE) 
rowlabT1<-C1[,1]
collabT1<-C1[,2]
p1<-randperm(length(rowlabT1),10000)
pre2<-predict(model,cbind(lncRNA_match_scale[rowlabT1[p1],],GExp_match_scale[collabT1[p1],]),probability=TRUE)
predictedY2<-attr(pre2,"probabilities")[,1]

CN<-which(abs(corMat)>=0 & abs(corMat)<0.001,arr.ind = TRUE) 
rowlabN<-CN[,1]
collabN<-CN[,2]
pN<-randperm(length(rowlabN),10000)
preN<-predict(model,cbind(lncRNA_match_scale[rowlabN[pN],],GExp_match_scale[collabN[pN],]),probability=TRUE)
predictedYN<-attr(preN,"probabilities")[,1]
# preN<-predict(model,Xn,probability=TRUE)

par(mfcol=c(1,4))
boxplot(predictedYP,notch=TRUE,cex.axis=1.5,ylab="Prediction Value",xlab=c("(0.7,1)"),cex.lab=1.5,col="gray")
#stripchart(predictedYP, vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg ="bisque",add = TRUE,cex=1.5)

boxplot(predictedY1,notch=TRUE,cex.axis=1.5,xlab=c("(0.5,0.7)"),cex.lab=1.5,col="gray")
#stripchart(predictedY1, vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg ="bisque",add = TRUE,cex=1.5)

boxplot(predictedY2,notch=TRUE,cex.axis=1.5,xlab=c("(0.001,0.1)"),cex.lab=1.5,col="gray")
#stripchart(predictedY2, vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg ="bisque",add = TRUE,cex=1.5)

boxplot(predictedYN,notch=TRUE,cex.axis=1.5,xlab=c("(0.0,0.00002)"),cex.lab=1.5,col="gray")
#stripchart(predictedYN, vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg ="bisque",add = TRUE,cex=1.5)


#lncRNA2Target low-throughput data test
load("lncRNA2Target-lowthroughput-data.RData")
overlap_lncRNA<-intersect(rownames(corMat),names(human_lncRNA))
overlap_lncRNA_target<-lapply(1:length(overlap_lncRNA),function(x) subset(human_lncRNA_data,lncRNA_name_from_paper==overlap_lncRNA[x])$Target_symbol_from_paper)
labR<-c(which(rownames(corMat) %in% "NEAT1"),which(rownames(corMat) %in% "HOTAIRM1"),which(rownames(corMat) %in% "HOTAIRM1"),which(rownames(corMat) %in% "HOTAIRM1"))
labC<-c(which(colnames(corMat) %in% "ADARB2"),which(colnames(corMat) %in% "HOXA1"),which(colnames(corMat) %in% "HOXA4"),which(colnames(corMat) %in% "EMX2"))

preL<-predict(model,cbind(lncRNA_match_scale[labR,],GExp_match_scale[labC,]),probability=TRUE)
predictedYL<-attr(preL,"probabilities")[,1]

load("lncRNA2Target-Highthroughput-data.RData")
overlap_lncRNA<-intersect(rownames(corMat),names(human_lncRNA))
overlap_lncRNA_target<-subset(human_lncRNA_data,lncrna_symbol==overlap_lncRNA)$gene_symbol
labR<-rep(which(rownames(corMat) %in% "NEAT1"),14)
for (i in 1:16) {labC[i]<-which(colnames(corMat) %in% intersect(overlap_lncRNA_target,colnames(corMat))[i])}

preH<-predict(model,cbind(lncRNA_ldata_BRCA[labR,],GExp_match[labC[-c(3,10)],]),probability=TRUE)
predictedYH<-attr(preH,"probabilities")[,1]
corMat_H<-diag(corMat[rowlabH,collabH])

# predicted the entile relationship matrix with row as lncRNAs, column as genes
corMat<-corMat[unique(rownames(corMat)),]
lncRNA_match_scale<-lncRNA_match_scale[unique(rownames(corMat)),]
GExp_match_scale[is.na(GExp_match_scale)]<-0
lncRNA_match_scale[is.na(lncRNA_match_scale)]<-0

nG<-ncol(corMat) # number of genes
nL<-nrow(corMat) # number of lncRNAs
len<-nL*nG
colabT<-vector();
for (i in 1:nG){colabT<-c(colabT,rep(i,nL))}
rowlabT<-rep(1:nL,nG)
t<-floor(len/10000)
prdY<-vector()
for (k in 1:10000) {
    ifelse (k==10000, {
      tk<-t*(k-1)+1;
      rowlab<-rowlabT[tk:len];collab<-colabT[tk:len];},{
        tk1<-t*(k-1)+1;tk2<-t*k;
        rowlab<-rowlabT[tk1:tk2];collab<-colabT[tk1:tk2];})       
    pre<-predict(model, cbind(lncRNA_ldata_BRCA[rowlab,],GExp_match[collab,]), probability=TRUE)
    predictedY<-attr(pre,"probabilities")[,1]
    prdY<-c(prdY,cbind(predictedY))
    cat(k,"\n")
}

preCOR<-matrix(prdY,nL,nG)
