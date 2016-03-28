library(stringr)
getAllExonsCoordiante<-function(){

  bed<-read.table("/home/limeng/splicingSNP/annotation/refgene.bed",header=FALSE,as.is=TRUE);

  chr<-bed[,1]
  exonsNumber<-bed[,10]
  totalExonsNumber<-sum(sapply(exonsNumber,function(x){ifelse(x>2,x,0)})  );
  
  resultExons<-vector(mode = "character", length = totalExonsNumber);
  
  transcriptBegPos<-bed[,2]
  transciptEndPos<-bed[,3]
  
  exonsLength<-bed[,11]
  exonsStart<-bed[,12]
  
  exonsLengthList<-strsplit(exonsLength,",");
  exonsStartList<-strsplit(exonsStart,",");
  j=1;
  
  for(i in 1:nrow(bed)){
    exonsBegPos<-transcriptBegPos[i]+as.numeric(exonsStartList[[i]])+1;
    exonsEndPos<-exonsBegPos+as.numeric(exonsLengthList[[i]]);
    
    transcriptExons<-str_c(chr[i],":",exonsBegPos,"-",exonsEndPos);
    
    #exclude the exons number <3
    if(length(transcriptExons)<=2)
      next;
    
    #exclude the first and last exons
    transcriptExons<-transcriptExons[-1];#exclude the first exon
    transcriptExons<-transcriptExons[-exonsNumber[i]]#exclude the last exon
    
    resultExons[j:(j+length(transcriptExons)-1 ) ]<-transcriptExons;
    
    j=j+length(transcriptExons)
      
  }
  
  cat(str_c("number of exons:",j))
  
  uniqueExons<-unique(resultExons);
  cat(str_c("number of unique exons:",length(uniqueExons),"\n"))
  uniqueExons<-uniqueExons[-length(uniqueExons)];
  
  
  #only keep chr[1-22,X,Y]
  uniqueExons1=uniqueExons[grepl("chr[1-9,x,y]:",uniqueExons,ignore.case = TRUE)]
  uniqueExons2=uniqueExons[grepl("chr1[0-9]:",uniqueExons,ignore.case = TRUE)]
  uniqueExons3=uniqueExons[grepl("chr2[0-2]:",uniqueExons,ignore.case = TRUE)]
  

  cat(str_c("number of unique exons only char1-22,x,y:",length(uniqueExons),"\n" ))
  cat(uniqueExons,file="allExons.csv",sep="\n")

}

setwd("/home/limeng/splicingSNP/1000genome/");
getAllExonsCoordiante();

