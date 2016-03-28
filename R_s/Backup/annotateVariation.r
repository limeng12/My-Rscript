library(stringr)

data<-read.table("/home/limeng/drugresistant/RNA-seq/mRNA/annotation/hg19refgene.gtf",header=FALSE,as.is=TRUE)

oldGeneId<-data[1,9];
exonIndex<-1;
for(i in 2:nrow(data)){
  if(data[i,3]!="exon")
    next;
  #cat(i)
  newGeneId<-data[i,9];
  
  ifelse(oldGeneId==newGeneId,{data[i,9]<-str_c(data[i,9],"exon number=",exonIndex);exonIndex<-exonIndex+1},exonIndex<-1);
  
  oldGeneId<-newGeneId;
  
  
}
