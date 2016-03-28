
featuresNames<-colnames(allDatasetTrain);

for(i in 2:ncol(allDatasetTrain)){
  
  currentFeature<-data.frame(allDatasetTrain[,1],allDatasetTrain[,i]);
  colnames(currentFeature)<-c("label","feature")
  
  
  lineLabel<-rep(TRUE,nrow(currentFeature));
  
  mccSumkknn<-0;
  mccSumsvm<-0;
  mccSumrandomforest<-0;
  
  for(j in 1:10){
    tlineLabel<-lineLabel;
    
    onePart<-floor(nrow(currentFeature)/10);
    
    tlineLabel[(1:onePart)*j]<-FALSE;
    testPart<-currentFeature[!tlineLabel,]
    trainPart<-currentFeature[tlineLabel,]
    
    resultkknn<-kknn(label~.,trainPart,testPart,k=7 ,kernel="biweight");
    mcckknn<-calculatingMCC(testPart$label,resultkknn$fitted.values );
    mccSumkknn<-mccSumkknn+mcckknn
    
    modelsvm<-ksvm(label~.,trainPart,kernel = "besseldot",kpar=list(sigma=1, order=0.001, degree=3)) 
    mccsvm<-calculatingMCC(testPart$label,predict(modelsvm,testPart[,"feature",drop=FALSE]) );    
    cat(mccsvm)
    
    mccSumsvm<-mccSumsvm+mccsvm
    
    modelrandomforest<-randomForest(formula=label~., data=trainPart);
    mccrandomforest<-calculatingMCC(testPart$label,predict(modelrandomforest,testPart) );
    mccSumrandomforest<-mccSumrandomforest+mccrandomforest;
    
    
  }
  
  cat(str_c("mcc of kknn=",format(mccSumkknn/10,digits=2),
            " mcc of svm=",  format(mccSumsvm/10,digits=2),
            " random forest=",format(mccSumrandomforest/10,digits=2),
            " feature:",featuresNames[i],"\n"));
  
}


pdf("qqplot and KS test.pdf")
par(mfrow=c(2,2))

#ss features
for(i in 2:13){
  testP<-ks.test(hgmdasass[,i],normalasass[,i])
  qqplot(normalasass[,i],hgmdasass[,i],ylab="hgmd",xlab="1000 genome",main=str_c("ss_",i-1," ks pvalue:",format(testP$p.value,digits=2)))
  abline(0,1)
}


#asa features
for(i in 14:16){
  testP<-ks.test(hgmdasass[,i],normalasass[,i])
  
  qqplot(normalasass[,i],hgmdasass[,i],ylab="hgmd",xlab="1000 genome",main=str_c("asa_",i-13," ks pvalue:",format(testP$p.value,digits=2)))
  abline(0,1)
  
}

#disorder features
for(i in 3:14){
  testP<-ks.test(hgmddisorder[,i],normaldisorder[,i])
  
  qqplot(normaldisorder[,i],hgmddisorder[,i],ylab="hgmd",xlab="1000 genome",main=str_c("disorder_",i-2," ks pvalue:",format(testP$p.value,digits=2)))
  abline(0,1)
  
}


#pfam features
testP<-ks.test(hgmdpfam[,ncol(hgmdpfam)],normalpfam[,ncol(normalpfam)])

qqplot(normalpfam[,ncol(normalpfam)],hgmdpfam[,ncol(hgmdpfam)],ylab="hgmd",xlab="1000 genome",main=str_c("pfam"," ks pvalue:",format(testP$p.value,digits=2)))
abline(0,1)




#ptm features
testP<-ks.test(rowSums(hgmdptm[,3:ncol(hgmdptm)]),rowSums(normalptm[,3:ncol(normalptm)]) )

less5hgmd<-rowSums(hgmdptm[,3:ncol(hgmdptm)])
less5normal<-rowSums(normalptm[,3:ncol(normalptm)])

x<-rowSums(hgmdptm[,3:ncol(hgmdptm)])
y<-rowSums(normalptm[,3:ncol(normalptm)]);

#qplot(x=x,y=y,stat="qq",xlab="hgmd",ylab="1000 genome",main=str_c("ptm"," ks pvalue:",testP$p.value) )
qqplot(y,x,ylab="hgmd",xlab="1000 genome",main=str_c("ptm"," ks pvalue:",format(testP$p.value,digits=2)),xlim=c(0,35),ylim=c(0,35) )
#abline(x=y)
abline(0,1)


dev.off();









