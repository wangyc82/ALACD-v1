crossvalSVM<-function(X,Y,fold,C,G)
{
  # this procedure is to implement SVM standard algorithm in cross-validation way;
  # X is inputs
  # Y is outputs
  # fold is number of K for K-fold cross-validation
  # C is the cost of constraints violation in Standard SVM model
  # G is rbf kernel parameter 
  prdY<-vector()
  
  library(pracma)
  library(e1071)
  len<-nrow(X) # the number of input
  
  temp<-randperm(len)
  Xr<-X[temp,];Yr<-Y[temp];
  #Xr<-X;Yr<-Y;
  t<-floor(len/fold);
  for (k in 1:fold) {
    if (k==fold) {tk<-t*(k-1)+1;Xtst<-Xr[tk:len,];Ytst<-Yr[tk:len];Xtrn<-Xr[-(tk:len),];Ytrn<-Yr[-(tk:len)];} else {tk1<-t*(k-1)+1;tk2<-t*k;Xtst<-Xr[tk1:tk2,];Ytst<-Yr[tk1:tk2];Xtrn<-Xr[-(tk1:tk2),];Ytrn<-Yr[-(tk1:tk2)];}       
    model<-svm(Xtrn,Ytrn,type="C-classification",kernel="radial",cost=C,gamma=G,class.weights=c(a=1, b=1),probability=TRUE)
    pre<-predict(model, Xtst,probability=TRUE)
    predictedY<-attr(pre,"probabilities")[,1]
    prdY<-c(prdY,cbind(predictedY))
    cat(k,"\n")
  }
  return(list(target=Yr, predY=prdY))
}
# usage
# crossV<-crossvalSVM(X,Ylab,fold,C,G)
#target<-crossV[[1]];prdY<-crossV[[2]]
# perf<-getperf(target,prediction) #get the criteria for performance