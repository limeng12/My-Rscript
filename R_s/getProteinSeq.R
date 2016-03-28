library(Biostrings)
library(ShortRead)
library(seqinr)

snps<-readLines("/home/limeng/splicingSNP/hgmdRawData/wrong exons.txt");
snpStr<-strsplit(snps,"\\:");

snpFrame<-data.frame( chr=as.character(sapply(snpStr,"[",1)),
                      pos=as.numeric(sapply(snpStr,"[",2)) ,
                      stringsAsFactors =FALSE);

mBed=new("BED");
mBed<-readFromBedFile(mBed,"/home/limeng/splicingSNP/annotation/refgene.bed");
id<-getGeneIdFromSnp(mBed,snpFrame);

nmMappingFile<-read.table(file="/home/limeng/splicingSNP/annotation/nm_np_mapping.txt",
                          header=FALSE,sep="\t",as.is=TRUE);

#npId<-nmMappingFile[which(nmMappingFile[,1] %in% id),2];
cat(id,file="nmId.csv",sep="\n");

x<-read.fasta("/home/limeng/splicingSNP/annotation/refgeneseqhg19.fa",
              as.string=TRUE,seqtype="AA");


#idNoVersion<-sapply( strsplit(as.character(id(x) ),"\\." ),"[",1 );
#result<-x[which(idNoVersion %in% npId) ];
#npProtein<-data.frame( np=sapply( strsplit(as.character(id(x)),"\\."),"[",1 ),
#                       seq=sread(x),stringsAsFactors = FALSE );

npProtein<-data.frame(
  np=sapply( strsplit(names(sapply(x,"[",1)),"\\.") ,"[",1), 
                      seq=unname(sapply(x,"[",1)),stringsAsFactors=FALSE );


nmProtein<-matrix(ncol=2,nrow=0);
cat("",file="result.txt",append=FALSE);  

id<-unique(id);

setwd(yaoNeutralPath);

for(i in 1:length(id) ) {
  
  npId<-nmMappingFile[ nmMappingFile[,1]==id[i] , 2 ];
  
  #nmId<-nmMappingFile[ nmMappingFile[,2]==npProtein[i,1] , 1 ];
  proteinSeq<-as.character(npProtein[ npProtein[,1]==npId,2 ]);

    
  if(length(npId)>1 ){
    cat("np ID length>1\n");
  }
  
  if(length(proteinSeq)==0|| (nchar(proteinSeq)<3) || is.na(proteinSeq) )
    next;
  
  for(j in 1:1) {
    
    #cat(">",file="result.txt",append=TRUE);  
    
    #cat(id[i],file="result.txt",append=TRUE); 
    #cat("\n",file="result.txt",append=TRUE);  
    
    #cat(npId[j],file="result.txt",append=TRUE); 
    #cat("\t",file="result.txt",append=TRUE);
    
    cat(proteinSeq,file=id[i],append=FALSE);  
    #cat("\n",file="result.txt",append=TRUE);  
    
  }
  
}

#write.table(nmProtein,file="result.txt");

