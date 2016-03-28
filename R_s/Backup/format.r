library(BSgenome.Hsapiens.UCSC.hg19)
setwd("/home/limeng/introSNPdisease")

datanormal<-read.table("normalsplicingsnpsf.msm",sep="\t",header=FALSE,as.is=TRUE);

snpSplit<-strsplit(datanormal[,1],":");
chr<-sapply(snpSplit,"[",1);
pos<-as.numeric(sapply(snpSplit,"[",2));

seq<-as.character(getSeq(Hsapiens,chr,pos,pos));
result<-cbind(datanormal[,1],seq,datanormal[,2],datanormal[,3],datanormal[,4]);
write.table(result,"nromalformat",row.names=FALSE,col.names=FALSE,quote=FALSE);
datahgmd<-read.table("hgmdsplicingsnpsf.msm",sep="\t",header=FALSE,as.is=TRUE);

snpSplit<-strsplit(datahgmd[,1],":")
chr<-sapply(snpSplit,"[",1)
pos<-as.numeric(sapply(snpSplit,"[",2))

seq<-as.character(getSeq(Hsapiens,chr,pos,pos));
result<-cbind(datahgmd[,1],seq,datahgmd[,2],datahgmd[,3],datahgmd[,4]);

write.table(result,"hgmdformat",row.names=FALSE,col.names=FALSE,quote=FALSE);



