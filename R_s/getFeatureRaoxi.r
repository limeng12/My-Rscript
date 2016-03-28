library(stringr);
library('BSgenome.Hsapiens.UCSC.hg19');

setwd("/home/limeng/Projects/matt");

#raoxiInput<-read.table("sigCoxTrain.txt",header=TRUE,as.is=TRUE,check.names=FALSE);
#names<-colnames(raoxiInput);
data<-read.table("exons.csv",header=FALSE,as.is=TRUE);
names<-data[,1];

centerExons<-sapply(strsplit(names[2:length(names)],"@"),"[",2);
exons<-strsplit(centerExons,":");


cat(centerExons,sep="\n",file="centerExons")

chr<-sapply(exons,"[",1);
beg<-as.numeric(sapply(exons,"[",2) );
end<-as.numeric(sapply(exons,"[",3) );
exonCenterPos<-round( (beg+end)/2);
strand<-sapply(exons,"[",4);

refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenterPos,exonCenterPos));

raoxi<-str_c(chr,":",exonCenterPos,"\t",refSeq,"\t","T","\t",strand,"\t",beg,"\t",end);

cat(raoxi,file="raoxi",sep="\n");

raoxiMap<-data.frame( id=str_c(chr,":",exonCenterPos),exon=names[2:length(names)] );

raoxiMap<-raoxiMap[!duplicated(raoxiMap[,1]),];
rownames(raoxiMap )<-raoxiMap[,1];


#setwd("/home/limeng/splicingSNP/code/raoxiFeature");
setwd("/home/limeng/splicingSNP/code/mattFeatures");

raoxiasass<-read.table("matt.asa.ss.max_prob.data",sep="\t",header=FALSE,as.is=TRUE);
raoxiasass<-raoxiasass[!duplicated(raoxiasass[,1]),];rownames(raoxiasass)<-raoxiasass[,1]
raoxidisorder<-read.table("matt.snp.ref.mut.seq.disorder-matrix",sep="\t",header=FALSE,as.is=TRUE);
raoxidisorder<-raoxidisorder[!duplicated(raoxidisorder[,1]),];rownames(raoxidisorder)<-raoxidisorder[,1];
raoxipfam<-read.table("matt.snp.ref.mut.seq.pfam-matrix",sep="\t",header=TRUE,as.is=TRUE);
raoxipfam<-raoxipfam[!duplicated(raoxipfam[,1]),];rownames(raoxipfam)<-raoxipfam[,1]
raoxiptm<-read.table("matt.snp.ref.mut.seq.ptm-matrix.dbPTM",sep="\t",header=TRUE,as.is=TRUE);
raoxiptm<-raoxiptm[!duplicated(raoxiptm[,1]),];rownames(raoxiptm)<-raoxiptm[,1];

raoxiPhylop<-read.table("matt.snp.ref.mut.seq.mutBed.phylop",sep="\t",header=FALSE);
raoxiPhylop<-raoxiPhylop[!duplicated(raoxiPhylop[,1]),];rownames(raoxiPhylop)<-(raoxiPhylop[,1]);


raoxiSNP<-intersect(raoxiasass[,1],raoxidisorder[,1]);
raoxiSNP<-intersect(raoxiSNP,raoxipfam[,1]);
raoxiSNP<-intersect(raoxiSNP,raoxiptm[,1]);
raoxiSNP<-intersect(raoxiSNP,raoxiPhylop[,1]);


raoxiData<-cbind(raoxiPhylop[raoxiSNP,2],raoxiasass[raoxiSNP,2:16],
                 raoxidisorder[raoxiSNP,3:14],raoxipfam[raoxiSNP,ncol(raoxipfam)],
                 rowSums(raoxiptm[raoxiSNP,3:ncol(raoxiptm)]) );

colnames(raoxiData)<-c("phylop","ss_1","ss_2","ss_3","ss_4","ss_5",
                       "ss_6","ss_7","ss_8","ss_9","ss_10","ss_11","ss_12",
                      "asa_1","asa_2","asa_3",
                      "disorder_1","disorder_2","disorder_3","disorder_4","disorder_5",
                      "disorder_6","disorder_7","disorder_8","disorder_9","disorder_10",
                      "disorder_11","disorder_12",
                      "pfam",
                      "ptm"
                      );

selectFeatureNames<-colnames(raoxiData);
#hgmdData<-hgmdData[hgmdData[,"phylop"]>hgmdPhylopThreshold,-1];
numberOfraoxiExons<-nrow(raoxiData);

#modelrandomforest<-randomForest(formula=label~., 
#                                data=alldataset[,c("label",selectFeatureNames)],
#                                ntree=500,proximity=TRUE,mtry=40,
#                                replace=FALSE,nodesize=20,maxnodes=70);

predictResult<-predict(modelrandomforestAll,
                       raoxiData[,selectFeatureNames],type="prob")[,1];

raoxiData<-cbind.data.frame(raoxiData,predictResult);
raoxiData<-cbind.data.frame(raoxiData,raoxiMap[rownames(raoxiData),2]);
colnames(raoxiData)[ncol(raoxiData)]<-"exonRegion";

#importantFeatures<-importance(modelrandomforest);
importantFeatures2<-importantFeatures[order(importantFeatures[,"MeanDecreaseGini"],decreasing=TRUE),];

selectNames<-names(importantFeatures2)[1:5];

setwd("/home/limeng/Projects/matt");
#setwd("/home/limeng/Projects/raoxi");

write.table(raoxiData[,c("exonRegion","predictResult",selectNames)],
          file="raoxiPredictResult.csv",row.names=FALSE,quote=FALSE,sep=",",col.names=c("exonRegion","predictResult",abv(selectNames) ) );


