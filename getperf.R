getperf<-function(target_class,prediction)
{
  #this program is to return the evalation criteria for target_class and the prediction
  library("ROCR")
  # target_class must be two level character vector
  cat("calculating the values for ROC curve \n")
  pred<- prediction(prediction, target_class)
  
  tpr <- performance(pred,"tpr");tpr<-unlist(slot(tpr, "y.values"));
  fpr <- performance(pred,"fpr");fpr<-unlist(slot(fpr, "y.values"));
  pre <- performance(pred,"prec");pre<-unlist(slot(pre, "y.values"));
  cat("calculating AUC \n") 
  auc <- performance(pred,"auc")
  auc<-unlist(slot(auc, "y.values"))
  
  cat("calculating AUPR \n")
  r<-vector()
  s<-vector()
  for (i in 1:length(tpr)-1) {r[i]<-(tpr[i+1]-tpr[i])*pre[i];s[i]<-0.5*(tpr[i+1]-tpr[i])*(pre[i+1]-pre[i]);}
  AUPR<-sum(r,na.rm=TRUE)+sum(s,na.rm=TRUE)
  
  P<-length(which(target_class %in% levels(as.factor(target_class))[1]));
  N<-length(which(target_class  %in% levels(as.factor(target_class))[1]));
  Sn<-1-fpr;
  Acc=(tpr*P+Sn*N)/(P+N);
  Fmeasure<-2*(tpr*Sn)/(tpr+Sn);
  cat("calculating the Acc (Accuracy), Sn (sensitivity), Sp (Specificity), Pre (precision), when F-measure reaches its maximun \n")
  y<-matrix(0,1,7)
  colnames(y)<-c("AUC","Acc","Sn","Sp","Pre","F-measure","AUPR")
  y[1]<-auc;y[2]<-Acc[which.max(Fmeasure)];y[3]<-tpr[which.max(Fmeasure)];y[4]<-Sn[which.max(Fmeasure)];
  y[5]<-pre[which.max(Fmeasure)];y[6]<-max(Fmeasure);y[7]<-AUPR
  
  getperf<-y;
}
