
pcaResult<-prcomp(alldatasetTrain[,-1]);
alldatasetPrcomp<-pcaResult$x;
alldatasetPrcomp<-cbind.data.frame(alldatasetTrain[,1],alldatasetPrcomp[,1:5]);

colnames(alldatasetPrcomp)<-c("label",str_c("prcomp",1:(ncol(alldatasetPrcomp) -1) ) );


sepLineNum<-round( nrow(alldatasetPrcomp)/3*2 );
dataTrainPrcomp<-alldatasetPrcomp[1:sepLineNum,];

dataTestPrcomp<-alldatasetPrcomp[( (sepLineNum+1):nrow(alldatasetTrain) ),];

modelrandomforest<-randomForest(formula=label~.,
                                data=dataTrainPrcomp,
                                ntree=100,mtry=40,proximity=TRUE,
                                replace=FALSE,nodesize=20,maxnodes=70);

#allDatasetTestROC<-
troc<-roc(dataTestPrcomp$label,
          predict(modelrandomforest,dataTestPrcomp,type="prob")[,1] );


par(mfrow=c(1,2));
plot(pcaResult$sdev/(sum(pcaResult$sdev)) ,
     ylab="variance contribution",
     xlab="feature index",
     main="PCA(svd) analysis" );



plot(troc,axes=FALSE,xlab="False Positive Rate",ylab="True Positive Rate");
text(0.5,1,labels=paste("Test on independent test set\n AUC=",
                        format(troc$auc,digits=4),"\n",sep=" ") );

p<-ggplot(alldatasetPrcomp)+geom_point(aes(x=prcomp1,y=prcomp2,color=label) );
print(p);

jpeg("pc1.jpg")
#plot(x=alldatasetPrcomp[alldatasetPrcomp[,1]=="NEUTRAL",2],y=rep(1,sum(alldatasetPrcomp[,1]=="NEUTRAL") ) );
p<-ggplot(alldatasetPrcomp)+geom_density(aes(x=prcomp1,color=label ) );

print(p);
dev.off();

jpeg("pc1Roc.jpg");
par(mfrow=c(1,1),mar=c(1,1,1,1));
troc<-roc(dataTestPrcomp$label,dataTestPrcomp$prcomp1   );

plot(troc,axes=FALSE,xlab="False Positive Rate",ylab="True Positive Rate");
text(0.5,1,labels=paste("AUC=",
                        format(troc$auc,digits=4),"\n",sep=" ") );

dev.off();


#dev.off();

