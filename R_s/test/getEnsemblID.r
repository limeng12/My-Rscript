setwd("/home/limeng/Projects/raoxi/raoxiEnsembl");
data<-read.table("seCoxSel_anno.txt",as.is=TRUE);

refseqIds<-data[,2];

refseqIdList<- str_split(refseqIds,","); 
refseqIdVec<- sapply(refseqIdList,"[",1); 

refseqIdVec<-refseqIdVec[!is.na(refseqIdVec)];


cat(refseqIdVec,file="miso_refseqIds",sep="\n"); 


data<-read.table("bioDBnet_db2db_150710124903_1187526809.txt",as.is=TRUE,sep="\t");

ensemblIdList<- str_split(data[,2],";"); 
ensemblIdVec<- sapply(ensemblIdList,"[",1); 
ensemblIdVec<-ensemblIdVec[!is.na(ensemblIdVec)];

cat(ensemblIdVec,file="miso_ensemblIds",sep="\n"); 


proteins<-readLines("raoxiEnsemblProteinSeq");

name<-"";
for(i in 1:length(proteins)){
  if(startsWith(proteins[i],">") ){
    name<-substr(proteins[i],2,nchar(proteins[i]) ); 
    
    next; 
  }
  
  cat(proteins[i],file=name,append=TRUE); 
}


for(i in 1:length(proteins)){
  if(startsWith(proteins[i],">") ){
    if(i!=1)
    cat("\n",file="proteinForPfam",append=TRUE)
    cat(proteins[i],file="proteinForPfam",append=TRUE); 
    cat("\n",file="proteinForPfam",append=TRUE)
    
    next; 
  }
  
  cat(proteins[i],file="proteinForPfam",append=TRUE); 
}

