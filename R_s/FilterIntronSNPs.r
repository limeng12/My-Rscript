#Not Used Anymore!!

#1 step 1 filter the duplicate snps.
#2 step 2 filter same exons
#3 step keep only GT-AG for positive strand

#snp structure #1chrosome 2.position 3.ref 4.alt 5.annotation

library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(plyr)

setwd("/home/limeng/splicingSNP/1000genomedata")

#keep the first snp;
removeDupSNPsIgnoreGenotype<-function(dataSNP){
  
  duplicateKeys<-duplicated(str_c(dataSNP[,1],dataSNP[,2])  )
  return(dataSNP[duplicateKeys,]);
}


#keep the first snp;
removeDupSNPs<-function(dataSNP){
  duplicateKeys<-duplicated(str_c(dataSNP[,1],dataSNP[,2],dataSNP[,3],dataSNP[,4])  )
  
  return(dataSNP[duplicateKeys,]);
  
}

a<-list()


#keep GT-AG CT-AC pattern
removeNonSplicingSNP<-function(dataSNP,befFileName){
  #pos<-as.numeric(sapply(chrPosList,"[",2) )
  #snp
  vcf<-with(dataSNP,GRanges(chr, IRanges(pos, pos),"*",ref,alt)  )
  
  bedData<-read.table(befFileName,sep="\t",header=FALSE,as.is=TRUE);
  colnames(bedData)<-c("chr","start","end","id","score","strand","thickStart","thickEnd","color","blockCount","blockSizes","blockStart");
  bed <- with(bedData, GRanges(chr, IRanges(start, end),strand,blockCount,blockSizes,blockStart));
  
  vcfhits<-findOverlaps(vcf,bed);
  #wrongAnnotationSNPs<-snp[-unique(queryHits(vcfhits))]
  #cat(str_c("!!wrong annotation snps:",wrongAnnotationSNPs,"\n") )
  
  vcfFrame<-as(vcf[queryHits(vcfhits)],"data.frame");
  bedFrame<-as(bed[subjectHits(vcfhits)],"data.frame");
  colnames(vcfFrame)<-str_c("snp",colnames(vcfFrame));
  colnames(bedFrame)<-str_c("bed",colnames(bedFrame));
  
  vcfbed<-cbind(vcfFrame,bedFrame);
  #snpchr
  
  #result<-apply(vcfbed,1,function(x){
    result<-apply(vcfbed,1,function(x){
    #result<-by(vcfbed,vcfbed[,c("snpseqnames","snpstart")],function(y){
    #cat("judge1\n")
    
      #return(apply(y,1,function(x){
        chr<-unique(x["snpseqnames"]);
        
        exonStartPos<-as.numeric(strsplit(x["bedblockStart"],",")[[1]] );  
        size<-as.numeric(strsplit(x["bedblockSizes"],",")[[1]]); 
        start<-as.numeric(strsplit(x["bedblockStart"],",")[[1]]);
        
        end<-start+size-1;  
        exonEndPos<-end;
        count<-as.numeric(x["bedblockCount"])
        strand<-x["bedstrand"];
        
        snppos<-as.numeric(x["snpstart"]);
        transcriptPos<-as.numeric(x["bedstart"]);
        transcriptEndPos<-as.numeric(x["bedend"]);
        relativeSnpPos<-snppos-transcriptPos-1;#conver to 0-based coordinate
        
        
        for(i in 1:count){
          disToStart<-relativeSnpPos-exonStartPos[i];
          disToEnd<-relativeSnpPos-exonEndPos[i];
          
          
          if((1==i)||(count==i)) #exclude the first and last exons
            next;
          
          
          if(  (abs(disToStart)>3)&&(abs(disToEnd)>3)  )
            next;
            
          if( (disToStart==-1)||(disToStart==-2) ){
              index=i;
            
          }else if( (disToEnd==1)||(disToEnd==2) ){
              index=i+1;
                
          }else
                next;
          
          
          if(strand=="+"){
            donor1<-exonEndPos[index-1]+transcriptPos+1+1;#G
            donor2<-exonEndPos[index-1]+transcriptPos+1+2;#T
            accept1<-exonStartPos[index]+transcriptPos+1-1;#G
            accept2<-exonStartPos[index]+transcriptPos+1-2;#A
            
            seq1<-as.character(getSeq(Hsapiens,rep(chr,4),c(donor1,donor2,accept1,accept2),c(donor1,donor2,accept1,accept2) )  );
            seq1<-toupper(seq1)
            
            if(   (seq1[1]!='G')||(seq1[2]!='T')||(seq1[3]!='G')||(seq1[4]!='A')   )
              next;
            
            return(x);
            
          }

          if(strand=="-"){
            accept1<-exonEndPos[index-1]+transcriptPos+1+1;#C
            accept2<-exonEndPos[index-1]+transcriptPos+1+2;#T
            donor1<-exonStartPos[index]+transcriptPos+1-1;#C
            donor2<-exonStartPos[index]+transcriptPos+1-2;#A
            
            
            seq2<-as.character(getSeq(Hsapiens,rep(chr,4),c(accept1,accept2,donor1,donor2),c(accept1,accept2,donor1,donor2) )  );
            seq2<-toupper(seq2)
            
            if(  (seq2[1]!='C')||(seq2[2]!='T')||(seq2[3]!='C')||(seq2[4]!='A')  )
              next;
            
            return(x);
            
         }
          
            
          
        }
      
        
      #}))
    
  })
  
  return(result)
}

setwd("/home/limeng/splicingSNP/1000genomedata")

data<-read.table("/home/limeng/splicingSNP/1000genomedata/1000genomesplcingIntronAF",as.is=TRUE,header=TRUE);
dataSnp<-data.frame(str_c("chr",data[,1]),data[,2],data[,4],data[,5]);
colnames(dataSnp)<-c("chr","pos","ref","alt")

bedPath<-"/home/limeng/splicingSNP/annotation/refgene.bed";

#dataSnp<-dataSnp[1:100,]
b<-removeNonSplicingSNP(dataSnp,bedPath)

data<-matrix(ncol=15,nrow=0);
for(i in 1:length(b)){
  if( length(b[[i]])!=0 )
    data<-rbind(data,b[[i]])
  
}

finalDataDub<-data[order(data[,1],data[,2]),];
write.csv(finalDataDub,file="with duplicate.csv",row.names=FALSE)

duplicatedKey<-duplicated(str_c(finalDataDub[,"snpseqnames"],finalDataDub[,"snpstart"],finalDataDub[,"snpref"],finalDataDub[,"snpalt"]) )

finalData<-finalDataDub[!duplicatedKey,]
write.csv(finalData,file="NeutralGTAGSNPAll.csv",quote=FALSE,row.names=FALSE)

#generate format for PSI
id<-str_c( finalData[,"snpseqnames"],":",as.numeric(finalData[,"snpstart"]),":",finalData[,"snpref"],">",finalData[,"snpalt"] )
psiFormat<-cbind( str_sub(finalData[,"snpseqnames"],4,nchar(finalData[,"snpseqnames"])),as.numeric(finalData[,"snpstart"]),id,finalData[,"snpref"],finalData[,"snpalt"]  );

write.table(psiFormat,file="ForPSI.csv",quote=FALSE,row.names=FALSE,sep="\t");




