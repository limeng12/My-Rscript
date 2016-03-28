library(stringr)
library(ggplot2)
library(kknn)
library(kernlab)
library(randomForest)
library(pROC)
library(gridBase)
library(grid)
library(gplots)
library(dplyr)
library(party)
#parameters

setwd("/Users/mengli/Documents/projects/splicingSNP/R/");
source("Global.r");
source("multiPlot.r");
#setwd("/home/limeng/splicingSNP/features/1000genome20");
setwd(workingDir);
clinvarExonIds<-readLines("clinvarExonId"); 
#plot(roc(allData$labelAll,allData$pfam));

options(warn=1);


loadData<-function(exclude_clinvar){
  
  setwd(HGMDFeaturePath);
  
  hgmdasass<-read.table("hgmddesp.asa.ss.max_prob.data",sep="\t",header=FALSE,as.is=TRUE);
  colnames(hgmdasass)<-c("snpId","ss_1","ss_2","ss_3","ss_4","ss_5",
                         "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                         "asa_1","asa_2","asa_3","exonId","asa_ave");
  
  hgmddisorder<-read.table("hgmddesp.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE,as.is=TRUE);
  colnames(hgmddisorder)<-c("snpId","region","disorder_1","disorder_2","disorder_3","disorder_4",
                            "disorder_5","disorder_6","disorder_7","disorder_8",
                            "disorder_9","disorder_10","disorder_11","disorder_12","exonId");
  proteinStarts<-as.numeric(sapply(strsplit(hgmddisorder[,"region"],":"),"[",2) );
  proteinEnds<-as.numeric(sapply(strsplit(hgmddisorder[,"region"],":"),"[",3) );
  proteinLength<-proteinEnds-proteinStarts+1
  hgmddisorder<-cbind(hgmddisorder,proteinLength);
  
  hgmdpfam<-read.table("hgmddesp.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE,as.is=TRUE);
  hgmdpfam<-data.frame(id=hgmdpfam[,"id"],pfam1=hgmdpfam[,"overlap"],
                       pfam2=hgmdpfam[,"coverage"],
                       exonId=hgmdpfam[,"exonId"],stringsAsFactors = FALSE ); 
  
  hgmdptm<-read.table("hgmddesp.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE,as.is=TRUE,quote=NULL);
  hgmdptm<-data.frame(id=hgmdptm[,"id"],ptm=rowSums(hgmdptm[,3:(ncol(hgmdptm)-1) ] ),
                      exonId=hgmdptm[,"exonId"] ,stringsAsFactors = FALSE); 
  
  hgmdPhylop<-read.table("hgmddesp.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE,as.is=TRUE);
  colnames(hgmdPhylop)<-c("snpId","phylop","r1","r2","r3","r4","exonId");
  
  hgmdData<-inner_join(hgmdasass,hgmddisorder,by=c("snpId"="snpId","exonId"="exonId") );
  hgmdData<-inner_join(hgmdData,hgmdpfam,by=c("snpId"="id","exonId"="exonId") );
  hgmdData<-inner_join(hgmdData,hgmdptm,by=c("snpId"="id","exonId"="exonId") );
  hgmdData<-inner_join(hgmdData,hgmdPhylop,by=c("snpId"="snpId","exonId"="exonId") );
  
  hgmdData<-hgmdData[!duplicated(hgmdData[,"snpId"]),];
  rownames(hgmdData)<-str_c(hgmdData[,"snpId"],rep("$",nrow(hgmdData)),
                            hgmdData[,"exonId"],rep("$",nrow(hgmdData)),
                            hgmdData[,"region"],rep("$HGMD",nrow(hgmdData)));
  if(exclude_clinvar){
    inClinvar<-is.element(hgmdData[,"exonId"],clinvarExonIds); 
    hgmdData<-hgmdData[!inClinvar,];
  }
  
  cat(hgmdData[,"exonId"],file=paste0(workingDir,"training_exonId.txt") ,sep="\n"); 
    
  allFeatureNames<-c("phylop",
                        "ss_1","ss_2","ss_3","ss_4","ss_5",
                        "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                        "asa_1","asa_2","asa_3",
                        "disorder_1","disorder_2","disorder_3","disorder_4",
                        "disorder_5","disorder_6","disorder_7","disorder_8",
                        "disorder_9","disorder_10","disorder_11","disorder_12",
                        "pfam1","pfam2","ptm","proteinLength");
  
  hgmdData<-hgmdData[,allFeatureNames];
  
  #hgmdData<-hgmdData[hgmdData[,"phylop"]>hgmdPhylopThreshold,-1];
  numberOfHGMDExons<-nrow(hgmdData);
  cat(paste0("number of HGMD exons:",numberOfHGMDExons,"\n") );
  
  setwd(yaoqiNeutralFeaturePath);#normalYaoqi
  
  normalasass<-read.table("normaldesp.asa.ss.max_prob.data",sep="\t",header=FALSE,as.is=TRUE);
  colnames(normalasass)<-c("snpId","ss_1","ss_2","ss_3","ss_4","ss_5",
                         "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                         "asa_1","asa_2","asa_3","exonId","asa_ave");
  
  normaldisorder<-read.table("normaldesp.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE,as.is=TRUE);
  colnames(normaldisorder)<-c("snpId","region","disorder_1","disorder_2","disorder_3","disorder_4",
                            "disorder_5","disorder_6","disorder_7","disorder_8",
                            "disorder_9","disorder_10","disorder_11","disorder_12","exonId");
  proteinStarts<-as.numeric(sapply(strsplit(normaldisorder[,"region"],":"),"[",2) );
  proteinEnds<-as.numeric(sapply(strsplit(normaldisorder[,"region"],":"),"[",3) );
  proteinLength<-proteinEnds-proteinStarts+1
  normaldisorder<-cbind(normaldisorder,proteinLength);
  
  
  normalpfam<-read.table("normaldesp.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE,as.is=TRUE);
  normalpfam<-data.frame(id=normalpfam[,"id"],pfam1=normalpfam[,"overlap"],
                         pfam2=normalpfam[,"coverage"],
                         exonId=normalpfam[,"exonId"],stringsAsFactors = FALSE ); 
  
  normalptm<-read.table("normaldesp.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE,as.is=TRUE,quote=NULL);
  normalptm<-data.frame(id=normalptm[,"id"],ptm=rowSums(normalptm[,3:(ncol(normalptm)-1) ] ),
                        exonId=normalptm[,"exonId"] ,stringsAsFactors = FALSE); 
  
  normalPhylop<-read.table("normaldesp.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE,as.is=TRUE);
  colnames(normalPhylop)<-c("snpId","phylop","r1","r2","r3","r4","exonId");
  
  normalData<-inner_join(normalasass,normaldisorder,by=c("snpId"="snpId","exonId"="exonId") );
  normalData<-inner_join(normalData,normalpfam,by=c("snpId"="id","exonId"="exonId") );
  normalData<-inner_join(normalData,normalptm,by=c("snpId"="id","exonId"="exonId") );
  normalData<-inner_join(normalData,normalPhylop,by=c("snpId"="snpId","exonId"="exonId") );
  
  normalData<-normalData[!duplicated(normalData[,"snpId"]),];
  rownames(normalData)<-str_c(normalData[,"snpId"],rep("$",nrow(normalData)),
                              normalData[,"exonId"],rep("$",nrow(normalData))
                              ,normalData[,"region"],rep("$NORMAL",nrow(normalData)));
  if(exclude_clinvar){
   inClinvar<-is.element(normalData[,"exonId"],clinvarExonIds);
   normalData<-normalData[!inClinvar,]; 
  }
  
  cat(normalData[,"exonId"],file=paste0(workingDir,"training_exonId.txt"),sep="\n",append=TRUE); 
  
  allFeatureNames<-c("phylop",
                     "ss_1","ss_2","ss_3","ss_4","ss_5",
                     "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                     "asa_1","asa_2","asa_3",
                     "disorder_1","disorder_2","disorder_3","disorder_4",
                     "disorder_5","disorder_6","disorder_7","disorder_8",
                     "disorder_9","disorder_10","disorder_11","disorder_12",
                     "pfam1","pfam2","ptm","proteinLength");
  
  normalData<-normalData[,allFeatureNames];
  ##filter by phylop score
  #normalData<-normalData[normalData[,"phylop"]<neutralPhylopThreshold,-1];
  numberOfNormalExons<-nrow(normalData);
  cat(paste0("number of neutral exons:",numberOfNormalExons,"\n") );
  
  allData<-rbind(normalData,hgmdData);
  labelAll<-as.factor( c(rep(0,nrow(normalData)),rep(1,nrow(hgmdData))  )  );
  allData<-cbind(labelAll,allData);
  
  #select mininum exons in each group
  if(nrow(normalData)<nrow(hgmdData) ){
    #hgmdData<-hgmdData[sample(1:nrow(hgmdData),nrow(normalData) ),];
    hgmdData<-hgmdData[(1:nrow(normalData) ),];
  }else{
    #normalData<-normalData[sample(1:nrow(hgmdData),nrow(hgmdData) ),];
    hgmdData<-hgmdData[(1:nrow(normalData) ),];
  }
  
  #0 represent negative set(neutral) and 1 represent positive set(HGMD).
  alldataset<-rbind(normalData,hgmdData)
  label<-as.factor(c(rep("NEUTRAL",nrow(normalData)),rep("HGMD",nrow(hgmdData) ) )  );
  alldataset<-cbind(label,alldataset )
  
  weka<-alldataset;
  #weka[,1]<-sapply(alldataset[,1],function(x){if("N"==x) return("NEUTRAL"); return("HGMD")  });
  
  write.csv(weka,"data.yaoqi.45vetebrate.phylop.csv",row.names=FALSE);
  
  #random shuffled
  set.seed(1);
  #set.seed(10); 
  alldataset<-alldataset[sample(1:nrow(alldataset)),];
  return(alldataset);

}

alldataset<-loadData(FALSE);
setwd(workingDir);

nLine<-as.integer(nrow(alldataset)/3*2);

#mission as NA not 0
alldatasetNA<-alldataset;
alldatasetNA[,"ptm"]<-sapply(alldataset[,"ptm"],
                             function(x){if(x==0) return(NA);return(x);});

alldatasetTrain<-alldataset[1:nLine,];
alldatasetTest<-alldataset[(nLine+1):nrow(alldataset),];

alldataset$label<-ordered(alldataset$label,levels=c("HGMD","NEUTRAL"));
alldatasetTrain$label<-ordered(alldatasetTrain$label,levels=c("HGMD","NEUTRAL"));
alldatasetTest$label<-ordered(alldatasetTest$label,levels=c("HGMD","NEUTRAL"));

write.table( alldataset,file=paste(workingDir,"alldataset.tsv",sep="" ),
             row.names=TRUE,col.names=T,quote=FALSE,sep="\t" );

#modelrandomforestAll<-randomForest(formula=label~.,
#                                data=alldataset,
#                                ntree=100,mtry=30,proximity=TRUE,
#                                replace=FALSE,nodesize=20,maxnodes=70);

modelrandomforestAll<-randomForest(formula=label~.,
                                data=alldataset,
                                ntree=100,mtry=12,proximity=TRUE,
                                replace=FALSE,nodesize=19); 

importantFeatures<-importance(modelrandomforestAll); 

save(alldataset,file="alldataset.Rdata"); 
#cmodelrf<-cforest(formula=label~.,alldataset,control = cforest_control(ntree = 50,mtry=12,replace=FALSE) ); 
#cimpFeatures<-varimp(cmodelrf,mincriterion=0,conditional=TRUE); 

featureNames<-dimnames(importantFeatures)[[1]][order(importantFeatures,decreasing=FALSE)]; 
gini<-importantFeatures[order(importantFeatures,decreasing=FALSE) ]; 


pdf("figure3_feature_importance.pdf"); 

importantFeatureNames<-featureNames[(length(featureNames)-9):length(featureNames)]; 
importantGini<-gini[(length(gini)-9):length(gini)]; 

par(mar=c(3,14,3,3),cex.lab=1,cex.axis=0.8); 
barplot(importantGini, horiz=TRUE, 
        names.arg=abv(importantFeatureNames), las=1, space=0.2,
        main="Features importance by Gini impuratity"
        ); 

dev.off(); 
cat(paste("current directory:",getwd(),"\n"))
save(modelrandomforestAll,file="model");

crossValidateRandomForest<-function(data,crossNum=10){
  allSampleIDs<-1:nrow(data);
  onePortionSize<-floor( length(allSampleIDs)/crossNum ); 
  
  sumRoc<-0;
  
  allProb<-c();
  allLabel<-c();
  
  for(i in 1:crossNum){
    set.seed(1000);
    testID<-sample(allSampleIDs,onePortionSize,replace=FALSE); 
    allSampleIDs<-setdiff(allSampleIDs,testID); 
    trainID<-setdiff(  1:nrow(data),  testID  ); 
    
    #modelrandomforest<-randomForest(formula=label~., 
    #                                data=data[trainID,],ntree=100,mtry=10, 
    #                                proximity=TRUE,replace=FALSE,nodesize=20,maxnodes=70 
    #);
    
    
    modelrandomforest<-randomForest(formula=label~., 
                                    data=data[trainID,],ntree=100,mtry=12,
                                    proximity=TRUE,replace=FALSE,nodesize=19
    );
    #troc<-(roc( (data[testID,])$label, predict(modelrandomforest,data[testID,],type="prob")[,1] )$auc) [1];
    prob<-predict(modelrandomforest,data[testID,],type="prob")[,1]; 
    label<-(data[testID,])$label; 
    
    allProb<-c(allProb,prob) ; 
    allLabel<-c(allLabel,as.character(label) ) ; 
    
    #sumRoc<-sumRoc+troc;
  } 
  return(data.frame(prob=allProb,label=allLabel  )  ); 
  #return (sumRoc/crossNum);
}


plotFunction<-function(predictValues,title){
  cat(paste("\ntrue positive for 0.9 threshold:") );
  cat(sum(predictValues[predictValues$label=="HGMD","prob"]>0.9)/sum(predictValues$label=="HGMD")  );
  
  cat(paste("\nfalse positive for 0.9 threshold:") );
  cat(sum(predictValues[predictValues$label=="NEUTRAL","prob"]>0.9)/sum(predictValues$label=="NEUTRAL")  );
  
  cat(paste("\ntrue positive for 0.8 threshold:") );
  cat(sum(predictValues[predictValues$label=="HGMD","prob"]>0.8)/sum(predictValues$label=="HGMD")  );
  
  cat(paste("\nfalse positive for 0.8 threshold:") );
  cat(sum(predictValues[predictValues$label=="NEUTRAL","prob"]>0.8)/sum(predictValues$label=="NEUTRAL")  );
  
  cat(paste("\ntrue positive for 0.5 threshold:") );
  cat(sum(predictValues[predictValues$label=="HGMD","prob"]>0.5)/sum(predictValues$label=="HGMD")  );
  
  cat(paste("\nfalse positive for 0.5 threshold:") );
  cat(sum(predictValues[predictValues$label=="NEUTRAL","prob"]>0.5)/sum(predictValues$label=="NEUTRAL")  );
  
  tprFun<-function(x){
    t<-sum(predictValues[predictValues$label=="HGMD","prob"]>x)/sum(predictValues$label=="HGMD");
    return(t);
  }
  
  fprFun<-function(x){
    t<-sum(predictValues[predictValues$label=="NEUTRAL","prob"]>x)/sum(predictValues$label=="NEUTRAL");
    return(t);
  }
  
  thresholds<-seq(-0.0001,1.0001,0.0001);
  tpr<-sapply(thresholds,tprFun);
  fpr<-sapply(thresholds,fprFun);
  
  cutoffs <- data.frame(cut=thresholds, tpr=tpr, fpr=fpr );
  write.csv(cutoffs,file="10_corss_validation_cutoff.csv",quote=FALSE,row.names=FALSE);
  
  #troc<-roc(allDatasetTestSelect$label,predictValues );
  #plot(troc);
  
  #aucValue<-troc$auc;
  aucValue<-mean(sample(predictValues[predictValues$label=="HGMD","prob"],1000000,replace=T) >
                   sample(predictValues[predictValues$label=="NEUTRAL","prob"],1000000,replace=T) );
  
  cat(paste("\nAUC=",aucValue,"\n")  );
  
  plot(x=fpr,y=tpr,xlab="False Positive Rate",ylab="True Positive Rate",type="l",main=title);
  abline(0,1);
  text(0.5,1,labels=paste("10 cross validation, AUC=", format(aucValue,digits=4),"\n",sep=" ") );
  
}

predictResultCV<-crossValidateRandomForest(alldataset,10);



#b<-ggplot(predictResult)+geom_density(aes(x=prob,color=label));
#densityNeutral<-density(predictResult[predictResult[,"label"]=="NEUTRAL",1],na.rm=TRUE );
#densityHGMD<-density(predictResult[predictResult[,"label"]=="HGMD",1],na.rm=TRUE );

#xMax<-max(c(densityNeutral$x,densityHGMD$x ) );xMin<-min(c(densityNeutral$x,densityHGMD$x ) ); 
#yMax<-max(c(densityNeutral$y,densityHGMD$y)  );yMin<-min(c(densityNeutral$y,densityHGMD$y)  ); 

#plot(densityNeutral$x,densityNeutral$y ,type="l",axes = F,
#     ylim=c(yMin,yMax),xlim=c(0,1),col="#E69F00",xlab="disease probability",ylab="density"); 

#lines(densityHGMD$x,densityHGMD$y,col="#000000"); 

#axis(1,cex.axis=0.7,mgp=c(0.5,0.5,0.5));
#axis(2,cex.axis=0.7,mgp=c(0.5,0.5,0.5));


bootstrapRandomForest<-function(data,bootstrapTime=100){
  allSampleIDs<-1:nrow(data);
  onePortionSize<-floor( length(allSampleIDs)/3 ); 
  
  sumRoc<-0;
  
  allProb<-c();
  allLabel<-c();
  
  for(i in 1:bootstrapTime){
    set.seed(1000);
    print(paste0("bootstrap number: ",i));
    testID<-sample(allSampleIDs,onePortionSize,replace=TRUE); 
    allSampleIDs<-setdiff(allSampleIDs,testID); 
    trainID<-setdiff(  1:nrow(data),  testID  ); 
    
    #modelrandomforest<-randomForest(formula=label~., 
    #                                data=data[trainID,],ntree=100,mtry=10,
    #                                proximity=TRUE,replace=FALSE,nodesize=20,maxnodes=70
    #);
    
    
    modelrandomforest<-randomForest(formula=label~., 
                                    data=data[trainID,],ntree=100,mtry=12,
                                    proximity=TRUE,replace=FALSE,nodesize=19
    );
    
    #troc<-(roc( (data[testID,])$label, predict(modelrandomforest,data[testID,],type="prob")[,1] )$auc) [1];
    prob<-predict(modelrandomforest,data[testID,],type="prob")[,1]; 
    label<-(data[testID,])$label; 
    
    allProb<-c(allProb,prob) ; 
    allLabel<-c(allLabel,as.character(label) ) ; 
    
    #sumRoc<-sumRoc+troc;
  } 
  return(data.frame(prob=allProb,label=allLabel  )  ); 
  #return (sumRoc/crossNum);
}

predictResultBoot<-predictResultCV<-crossValidateRandomForest(alldataset,100);

#pdf("figure9.pdf",width=8,height=5);
jpeg("supplement_figure2.jpeg",width=1200,height=600);
par(mfrow=c(1,2) );

plotFunction(predictResult,"10-cross validation");
plotFunction(predictResultBoot,"100 times bootstrap validation");

dev.off();
