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

testUsingPSI<-function(){
  
  setwd(psiNeutralFeaturePath);
  
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
  alldatasetPSI<-normalData;
  
  write.csv(alldatasetPSI[,1:2],file="psiData.csv",quote=FALSE);
  
  
  label<-rep("NEUTRAL",nrow(alldatasetPSI) );
  alldatasetPSI<-cbind(label,normalData);
 
  #random shuffled
  set.seed(1);
  alldatasetPSI<-alldatasetPSI[sample(1:nrow(alldatasetPSI)),];
  
  #nLine<-as.integer(nrow(alldatasetPSI)/3*2);
  
  target<-c("label","phylop","ss_1","ss_2","ss_3","ss_4",
            "ss_5","ss_6","ss_7","ss_8","ss_9",
            "ss_10","ss_11","ss_12",
            "asa_1","asa_2","asa_3",
            "disorder_1","disorder_2","disorder_3",
            "disorder_4","disorder_6","disorder_7","disorder_8",
            "disorder_9","disorder_10","disorder_11","disorder_12",
            "pfam");
  
  #allDatasetTrain<-alldatasetPSI[1:nLine,target];
  #allDatasetTest<-alldatasetPSI[(nLine+1):nrow(alldataset),target];
  
  
  #modelrandomforest<-randomForest(formula=label~., data=alldataset,
  #                                ntree=500,mtry=40,proximity=TRUE,
  #                                replace=FALSE,nodesize=20,maxnodes=70);
  
  predictPSI<<-predict(modelrandomforestAll,alldatasetPSI,type="prob")[,1];
  #hist(predictPSI);
 
  #troc<-roc(alldatasetPSI$label, predictPSI);
  
  #plot(troc,main=paste("auc=",format(troc$auc,digits=4),"method=random forest",sep=" ") );
  
  #first 1000 is neutral, last 1000 is HGMD
  idAFPSI<-read.csv(file=paste(psiFilesDir,"mappingFile.csv",sep=""),
                    header=FALSE,as.is=TRUE);
  
  rownames(idAFPSI)<-idAFPSI[,1];
  #b<-(idAFPSI[names(a)[1:2000],])
  
  #cat( names(a),file="predictNames.csv",sep="\n" )
  neutralIDPsi<-sapply(strsplit(rownames(alldatasetPSI[alldatasetPSI[,1]=="NEUTRAL",]) ,"\\$") ,"[",1);
  neutralID<-rownames(alldatasetPSI[alldatasetPSI[,1]=="NEUTRAL",]) ;
  
  #hgmdID<-rownames(alldatasetPSI[alldatasetPSI[,1]=="HGMD",]);
  neutralTable<-cbind(idAFPSI[neutralIDPsi,2:3],predictPSI[neutralID] );
  #hgmdTable<-cbind(predictPSI[hgmdID]);
  
  neutralTableValidation<-cbind(idAFPSI[neutralIDPsi,1:4],predictPSI[neutralID]  );
  write.csv(neutralTableValidation,file="neutralTableValidation.csv");
  
  colnames(neutralTable)<-c("AF","PSI","score");
  #colnames(hgmdTable)<-c("score");
  
  neutralTable<-neutralTable[!is.na(neutralTable[,1]),];
  neutralTable<-transform(neutralTable, AF = as.numeric(AF),
                          PSI = as.numeric(PSI),score = as.numeric(score)  );
  
  #hgmdTable<-transform(hgmdTable, score = as.numeric(score)  );
  
 pdf("using our model on 1000genome.pdf");
 neutralPsi5<-subset(neutralTable,abs(PSI)>5);
 
 p<-ggplot( neutralPsi5 )+geom_histogram(aes(x=score),binwidth =0.05)+
   ggtitle("using our model on 1000genome (|PSI|>5% MAF>5%)")+xlab("disease causing probability");
 
 print(p);
 #multiplot(p1,cols=1 );
 dev.off();
 
 
  write.csv(neutralTable,file="neutral.csv",row.names=FALSE,quote=FALSE); 
  
  psiBigger20<-neutralTable[abs(neutralTable[,"PSI"])>20,];
  psiBigger10Smaller20<-neutralTable[( abs(neutralTable[,"PSI"])<20 )&(abs(neutralTable[,"PSI"])>10  ),];
  psiSmaller10<-neutralTable[ abs(neutralTable[,"PSI"])<10,];
  
  cat(str_c("|PSI|>20 number of samples score>0.9=" ,
            nrow(psiBigger20[psiBigger20[,"score"]>0.9,] )/nrow(psiBigger20) ,"\n")  );
  
  cat(str_c("10<|PSI|<20 number of samples score>0.9=" , 
            nrow(psiBigger10Smaller20[psiBigger10Smaller20[,"score"]>0.9,] )/nrow(psiBigger10Smaller20) ,"\n")  );
  
  cat(str_c("|PSI|<10 number of samples score>0.9=" , 
            nrow(psiSmaller10[psiSmaller10[,"score"]>0.9,] )/nrow(psiSmaller10) ,"\n")  );
  

}

testUsingPSI();



