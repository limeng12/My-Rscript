library(stringr)
library(ggplot2)
library(kknn)
library(kernlab)
library(randomForest)
library(pROC)
library(gridBase)
library(grid)
library(gplots)
#parameters:
neutralPhylopThreshold<-(0.05);
#parameters:
hgmdPhylopThreshold<-(0.35);


setwd("/home/limeng/splicingSNP/R/");
source("multiPlot.r");
setwd("/home/limeng/splicingSNP/features/1000genome20");

#plot(roc(allData$labelAll,allData$pfam));


setwd("/home/limeng/splicingSNP/code")

hgmdasass<-read.table("hgmd.asa.ss.max_prob.data",sep="\t",header=FALSE,as.is=TRUE);
hgmdasass<-hgmdasass[!duplicated(hgmdasass[,1]),];rownames(hgmdasass)<-hgmdasass[,1]
hgmddisorder<-read.table("hgmd.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE,as.is=TRUE);
hgmddisorder<-hgmddisorder[!duplicated(hgmddisorder[,1]),];rownames(hgmddisorder)<-hgmddisorder[,1];
hgmdpfam<-read.table("hgmd.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE,as.is=TRUE);
hgmdpfam<-hgmdpfam[!duplicated(hgmdpfam[,1]),];rownames(hgmdpfam)<-hgmdpfam[,1]
hgmdptm<-read.table("hgmd.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE,as.is=TRUE);
hgmdptm<-hgmdptm[!duplicated(hgmdptm[,1]),];rownames(hgmdptm)<-hgmdptm[,1];

hgmdPhylop<-read.table("hgmd.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE);
hgmdPhylop<-hgmdPhylop[!duplicated(hgmdPhylop[,1]),];rownames(hgmdPhylop)<-(hgmdPhylop[,1]);


hgmdSNP<-intersect(hgmdasass[,1],hgmddisorder[,1]);
hgmdSNP<-intersect(hgmdSNP,hgmdpfam[,1]);
hgmdSNP<-intersect(hgmdSNP,hgmdptm[,1]);
hgmdSNP<-intersect(hgmdSNP,hgmdPhylop[,1]);


hgmdData<-cbind(hgmdPhylop[hgmdSNP,2],hgmdasass[hgmdSNP,2:16],
                hgmddisorder[hgmdSNP,3:14],hgmdpfam[hgmdSNP,ncol(hgmdpfam)],
                rowSums(hgmdptm[hgmdSNP,3:ncol(hgmdptm)]) )
colnames(hgmdData)<-c("phylop","ss_1","ss_2","ss_3","ss_4",
                      "ss_5","ss_6","ss_7","ss_8",
                      "ss_9","ss_10","ss_11","ss_12",
                      "asa_1","asa_2","asa_3",
                      "disorder_1","disorder_2","disorder_3",
                      "disorder_4","disorder_5","disorder_6",
                      "disorder_7","disorder_8","disorder_9",
                      "disorder_10","disorder_11","disorder_12",
                      "pfam",
                      "ptm"
)

#hgmdData<-hgmdData[hgmdData[,"phylop"]>hgmdPhylopThreshold,-1];
numberOfHGMDExons<-nrow(hgmdData)

setwd("/home/limeng/splicingSNP/code/normalAF5poPSI5/");
#setwd("/home/limeng/splicingSNP/code/normalAF5PSI5/");
#setwd("/home/limeng/splicingSNP/code/normalYaoqi/");#normalYaoqi


normalasass<-read.table("normal.asa.ss.max_prob.data",sep="\t",header=FALSE);
normalasass<-normalasass[!duplicated(normalasass[,1]),];rownames(normalasass)<-(normalasass[,1]);
normaldisorder<-read.table("normal.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE);
normaldisorder<-normaldisorder[!duplicated(normaldisorder[,1]),];rownames(normaldisorder)<-(normaldisorder[,1]);            
normalpfam<-read.table("normal.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE);
normalpfam<-normalpfam[!duplicated(normalpfam[,1]),];rownames(normalpfam)<-(normalpfam[,1]);
normalptm<-read.table("normal.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE);
normalptm<-normalptm[!duplicated(normalptm[,1]),];rownames(normalptm)<-(normalptm[,1]);

normalPhylop<-read.table("normal.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE);
normalPhylop<-normalPhylop[!duplicated(normalPhylop[,1]),];rownames(normalPhylop)<-(normalPhylop[,1]);


#setwd("/home/limeng/splicingSNP/features/AF5PSI5hgmd/")
setwd("/home/limeng/splicingSNP/features/AF5poPSI5hgmd/");
#setwd("/home/limeng/splicingSNP/features/YaoqiExonshgmd/");


normalSNP<-intersect(normalasass[,1],normaldisorder[,1])
normalSNP<-intersect(normalSNP,normalpfam[,1])
normalSNP<-intersect(normalSNP,normalptm[,1])
normalSNP<-intersect(normalSNP,normalPhylop[,1])


normalData<-cbind(normalPhylop[normalSNP,2],normalasass[normalSNP,2:16],
                  normaldisorder[normalSNP,3:14],normalpfam[normalSNP,ncol(normalpfam)],
                  rowSums(normalptm[normalSNP,3:ncol(normalptm)]) )

colnames(normalData)<-c("phylop","ss_1","ss_2","ss_3",
                        "ss_4","ss_5","ss_6","ss_7","ss_8",
                        "ss_9","ss_10","ss_11","ss_12",
                        "asa_1","asa_2","asa_3",
                        "disorder_1","disorder_2","disorder_3",
                        "disorder_4","disorder_5","disorder_6",
                        "disorder_7","disorder_8","disorder_9",
                        "disorder_10","disorder_11","disorder_12",
                        "pfam",
                        "ptm"
)


##filter by phylop score
#normalData<-normalData[normalData[,"phylop"]<neutralPhylopThreshold,-1];
numberOfNormalExons<-nrow(normalData);

#allData<-rbind(normalData,hgmdData);
#labelAll<-as.factor( c(rep(0,nrow(normalData)),rep(1,nrow(hgmdData))  )  );
#allData<-cbind(labelAll,allData);

if(nrow(normalData)<nrow(hgmdData) ){
  #hgmdData<-hgmdData[sample(1:nrow(hgmdData),nrow(normalData) ),];
  hgmdDataLast1000<-hgmdData[(nrow(hgmdData)-1000+1 ):nrow(hgmdData) ,];
  normalDataFirst1000<-normalData[1:1000,]
}else{
  #normalData<-normalData[sample(1:nrow(hgmdData),nrow(hgmdData) ),];
  hgmdData<-hgmdData[(1:nrow(normalData) ),];
}


#0 represent negative set(neutral) and 1 represent positive set(HGMD).
alldatasetPSI<-rbind(normalDataFirst1000,hgmdDataLast1000)
label<-as.factor(c(rep("NEUTRAL",nrow(normalDataFirst1000)),rep("HGMD",nrow(hgmdDataLast1000)  ) )  );
alldatasetPSI<-cbind(label,alldatasetPSI  )  

write.csv(alldatasetPSI[,1:2],file="psiData.csv",quote=FALSE);

#weka<-alldataset;
#weka[,1]<-sapply(alldataset[,1],function(x){if("N"==x) return("NEUTRAL"); return("HGMD")   });

#write.csv(weka,"data.yaoqi.45vetebrate.phylop.csv",row.names=FALSE);

#random shuffled
set.seed(1);
alldatasetPSI<-alldatasetPSI[sample(1:nrow(alldatasetPSI)),];

nLine<-as.integer(nrow(alldatasetPSI)/3*2);

target<-c("label","phylop","ss_1","ss_2","ss_3","ss_4","ss_5","ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
          "asa_1","asa_2","asa_3",
          "disorder_1","disorder_2","disorder_3","disorder_4","disorder_6","disorder_7","disorder_8","disorder_9","disorder_10","disorder_11","disorder_12",
          "pfam");

#allDatasetTrain<-alldatasetPSI[1:nLine,target];
#allDatasetTest<-alldatasetPSI[(nLine+1):nrow(alldataset),target];


modelrandomforest<-randomForest(formula=label~., data=alldataset,ntree=500,mtry=40,proximity=TRUE,replace=FALSE,nodesize=20,maxnodes=70);
predictPSI<-predict(modelrandomforest,alldatasetPSI,type="prob")[,1];
troc<-roc(alldatasetPSI$label, predictPSI);

plot(troc,main=paste("auc=",format(troc$auc,digits=4),"method=random forest",sep=" ") );

#first 1000 is neutral, last 1000 is HGMD
idAFPSI<-read.csv(file="/home/limeng/splicingSNP/1000genome/WholeExonIntron20/PSIResult/AF5population/mappingFile.csv",header=FALSE,as.is=TRUE);
rownames(idAFPSI)<-idAFPSI[,1];
#b<-(idAFPSI[names(a)[1:2000],])

#cat( names(a),file="predictNames.csv",sep="\n" )
neutralID<-rownames(alldatasetPSI[alldatasetPSI[,1]=="NEUTRAL",]);
hgmdID<-rownames(alldatasetPSI[alldatasetPSI[,1]=="HGMD",]);
neutralTable<-cbind(idAFPSI[neutralID,2:3],predictPSI[neutralID]  );
hgmdTable<-cbind(predictPSI[hgmdID]);

colnames(neutralTable)<-c("AF","PSI","score");
colnames(hgmdTable)<-c("score");

neutralTable<-neutralTable[!is.na(neutralTable[,1]),];

neutralTable<-transform(neutralTable, AF = as.numeric(AF),PSI = as.numeric(PSI),score = as.numeric(score)    );
hgmdTable<-transform(hgmdTable, score = as.numeric(score)  );


write.csv(neutralTable,file="neutral.csv",row.names=FALSE,quote=FALSE);
write.csv(hgmdTable,file="hgmd.csv",row.names=FALSE,quote=FALSE);

smoothScatter(abs(neutralTable[,2]),neutralTable[,3]);

neutralTableBigger0<-neutralTable[neutralTable[,3]<0.05,]
#neutralTableBigger0<-neutralTableBigger0[neutralTableBigger0[,1]>0.2,]


smoothScatter(abs(neutralTableBigger0[,2]),neutralTableBigger0[,3],xlab="|dPSI|",ylab="prediction probability",main="prediction score < 0.05");






