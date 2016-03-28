setwd("/home/limeng/splicingSNP/1000genome/WholeExonIntron20/evolution comparison/")

b=new("BED");
mBed<-readFromBedFile(b,"/home/limeng/splicingSNP/annotation/refgene.bed");
dataFileAF5<-read.table(file="AF5.rand.forPSI",as.is=TRUE);
snpAF5<-data.frame(chr=str_c("chr",dataFileAF5[,1]),pos=dataFileAF5[,2],stringsAsFactors =FALSE);

AF5Exons<-getNearestExons(mBed,snpAF5,3);AF5Exons<-unique(AF5Exons);
cat(AF5Exons,file="AF5Exons.csv",sep="\n");

dataFileAF10<-read.table("AF10.rand.forPSI",as.is=TRUE);
snpAF10<-data.frame(chr=str_c("chr",dataFileAF10[,1]),pos=dataFileAF10[,2],stringsAsFactors =FALSE);

AF10Exons<-getNearestExons(mBed,snpAF10,3);AF10Exons<-unique(AF10Exons);
cat(AF5Exons,file="AF10Exons.csv",sep="\n");


setwd("/home/limeng/splicingSNP/code");

dataAF5<-read.table("AF5Exons.csv.phylop",header=FALSE,as.is=TRUE);
dataAF5<-cbind(dataAF5[,1],rep("dataAF5",nrow(dataAF5)));
dataAF5<-unname(dataAF5);

dataAF10<-read.table("AF10Exons.csv.phylop",header=FALSE,as.is=TRUE);
dataAF10<-cbind(dataAF10[,1],rep("dataAF10",nrow(dataAF10)));
dataAF10<-unname(dataAF10);

allData<-rbind(dataAF5,dataAF10);
allData<-data.frame(evolution=as.numeric(allData[,1]),label=allData[,2]);
colnames(allData)<-c("evolution.score","label");

pd<-ggplot(allData);
pd<-pd+geom_density(aes(x=evolution.score,y=..density..,colour=label),position="identity" );

print(pd);

