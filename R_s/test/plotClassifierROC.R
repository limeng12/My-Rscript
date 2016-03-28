target<-c("label","phylop","ss_2","ss_3","ss_5","ss_7","ss_8",
          "asa_1",
          "disorder_2","disorder_3","disorder_12",
          "pfam","ptm");

target<-c("label","phylop","ss_1","ss_2","ss_3",
          "ss_4","ss_5","ss_6","ss_7","ss_8",
          "ss_9","ss_10","ss_11","ss_12",
          "asa_1","asa_2","asa_3",
          "disorder_1","disorder_2","disorder_3","disorder_4",
          "disorder_6","disorder_7","disorder_8","disorder_9",
          "disorder_10","disorder_11","disorder_12",
          "pfam","ptm");

allDatasetTrain<-alldataset[1:nLine,target];
allDatasetTest<-alldataset[(nLine+1):nrow(alldataset),target];


#modelsvm<-ksvm(label~.,allDatasetTrain,kpar=list(sigma=0.05),kernel = "rbfdot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,sep=""));
#cat(str_c("svm mcc rbfdot:",mccall));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "polydot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,"kernal=","polydot",sep=""));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "vanilladot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,"kernal=","vanilladot",sep=" "));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "tanhdot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,"kernal=","tanhdot",sep=""));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "laplacedot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,"kernal=","laplacedot",sep=""));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "besseldot") ;
#mccall<-calculatingMCC(allDatasetTest$label,predict(modelsvm,allDatasetTest) ); 
#cat(str_c("svm mcc besseldot:",mccall));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "splinedot") ;
#mccall<-calculatingMCC(allDatasetTest$label,predict(modelsvm,allDatasetTest) );
#cat(str_c("svm mcc splinedot:",mccall));

#modelsvm<-ksvm(label~.,allDatasetTrain,kernel = "stringdot") ;
#mccall<-calculatingMCC(allDatasetTest$label,predict(modelsvm,allDatasetTest) );
#cat(str_c("svm mcc stringdot:",mccall));
#pdf("ROC and AUC.pdf")

#modelsvm<-ksvm(label~.,allDatasetTrain,kpar=list(sigma=0.05),kernel = "rbfdot",prob.model = TRUE) ;
#troc<-roc(allDatasetTest$label,predict(modelsvm,allDatasetTest,type="probabilities")[,1]);
#plot(troc,main=paste("auc=",troc$auc,sep=""));
#modelrandomforest<-randomForest(formula=label~., data=allDatasetTrain,ntree=500,mtry=40,proximity=TRUE,replace=FALSE,nodesize=20,maxnodes=70);

modelrandomforest<-randomForest(formula=label~., 
                                data=allDatasetTrain,ntree=500,
                                mtry=40,proximity=TRUE,replace=FALSE,nodesize=20,maxnodes=70);
troc<-roc(allDatasetTest$label,
          predict(modelrandomforest,allDatasetTest,type="prob")[,1] );

plot(troc,main=paste("auc=",format(troc$auc,digits=4),"method=random forest",sep=" ") );

#dev.off();
library(ROCR)

pred <- prediction( predict(modelrandomforest,allDatasetTest,type="prob")[,1], allDatasetTest$label)
perf <- performance(pred,"tpr","fpr")

cutoffs <- data.frame(cut=format(perf@alpha.values[[1]],digits=3), 
                      fpr=format(perf@y.values[[1]],digits=3), tpr=format(perf@x.values[[1]],digits=3)  )

#cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@y.values[[1]], tpr=perf@x.values[[1]]  )

colnames(cutoffs)<-c("threshold","false positive rate","true positive rate");
write.csv(cutoffs,file="roc.csv",quote=FALSE,row.names=FALSE);

#plot(cutoffs[,2],cutoffs[,3])


