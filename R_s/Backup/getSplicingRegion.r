

generatingIntronSplicingRegion<-function(bedFileName,rightIntronLength,leftExonLength,leftIntronLength,rightExonLength){
  bedData<-read.table(bedFileName,header=FALSE,as.is=TRUE);
  colnames(bedData)<-c("chr","start","end","id","score","strand","thickStart","thickEnd","color","blockCount","blockSizes","blockStart");
  
  #intronSplicingRegion<-list();
  leftExonLength<-leftExonLength-1;
  rightExonLength<-rightExonLength-1;
  
  chrs<-bedData[,"chr"];
  exonStartPoss<-strsplit(bedData[,"blockStart"],",");  
  sizes<-strsplit(bedData[,"blockSizes"],","); 
  counts<-as.numeric(as.character(bedData[,"blockCount"]));
  transcriptPoss<-bedData[,"start"];
  
  #rowNumber<-sum(counts);
  rowNumber<-0;
  for(i in 1:length(counts)){
    count<-counts[i]
    if(count<3)
      next;
    
    rowNumber=rowNumber+counts[i]-2;
    
  }
  
  
  intronSplicingRegion<- matrix(ncol=3,nrow=rowNumber*2);
  index=1;
  
  for(i in 1:nrow(bedData) ){
    #x<-bedData[i,];
    if(i%%100==1)
      cat(str_c(i,"\n"));
    
    chr<-chrs[i];
    start<-as.numeric(exonStartPoss[[i]])
    size<-as.numeric(sizes[[i]])
    count<-counts[i]
    transcriptStartPos<-transcriptPoss[i]
        
    end<-start+size-1;  
    exonEndPos<-end;
    
    transcriptStart<-start+transcriptStartPos+1
    transcriptEnd<-transcriptStartPos+end+1;
    
    #cat(str_c("exon count:",count,"\n"))
    
    for(j in 1:count){
      
      trightIntronLength<-rightIntronLength
      tleftExonLength<-leftExonLength
      tleftIntronLength<-leftIntronLength
      trightExonLength<-rightExonLength
      
      if(1==j)
        next;
      if(count==j)
        next;
      lastIntronLength=start[j]-end[j-1]+1;
      thisExonsLength=size[j];
      nextIntronLength=start[j+1]-end[j]+1;
      
      if(tleftIntronLength>lastIntronLength )
        tleftIntronLength=lastIntronLength
      
      if(trightExonLength>thisExonsLength)
        trightExonLength=thisExonsLength
      
      if(tleftExonLength>thisExonsLength)
        tleftExonLength=thisExonsLength  
      
      if(trightIntronLength>nextIntronLength)
        trightIntronLength=nextIntronLength
      
      
      stopifnot(leftExonLength<=2);
      #here -1 means convert to bed format
      intronSplicingRegion[index,] <- c(chr,transcriptStart[j]-1-tleftIntronLength,transcriptStart[j]-1+trightExonLength);
      index=index+1;
      intronSplicingRegion[index,] <- c(chr,transcriptEnd[j]-1-tleftExonLength,transcriptEnd[j]-1+trightIntronLength);
      index=index+1;
      
      #intronSplicingRegion[index,] <- c(chr,transcriptStart[j]-tleftIntronLength,transcriptStart[j]+trightExonLength);
      #index=index+1;
      #intronSplicingRegion[index,] <- c(chr,transcriptEnd[j]-tleftExonLength,transcriptEnd[j]+trightIntronLength);
      #index=index+1;
      
    }
    
  }
  
  #return(intronSplicingRegion);
  
  #post sort and processing,sort and annotation
  intronBed<-data.frame(intronSplicingRegion[,1],as.numeric(intronSplicingRegion[,2]),
                        as.numeric(intronSplicingRegion[,3]),stringsAsFactors =FALSE);
  
  intronBedunique<-unique(intronBed);
  intronBedSort<-intronBedunique[order(intronBedunique[,1],intronBedunique[,2],intronBedunique[,3]),  ]
  
  #intronBedSort<-cbind(intronBedSort)
  #colnames(intronBedSort)<-c("chrom","chromStart","chromEnd","name","score",
  #"strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts");
  colnames(intronBedSort)<-c("chrom","chromStart","chromEnd")
  
  #delete "chr" in exons
  cat(str_c("number of bed's chrosome doesn't begin with 0",sum(!str_detect(intronBedSort[,1],"chr") ) )  );
    
  if(sum(!str_detect(intronBedSort[,1],"chr") )==0)
    intronBedSort[,1]<-str_sub(intronBedSort[,1],4)
  
  return(intronBedSort);
  
}


#intronBedSort<-generatingIntronSplicingRegion("refgene.bed",2,0);

#options("scipen"=1000, "digits"=4);
#write.table(intronBedSort,file="splicingRegionBed2.bed",row.names=FALSE,quote=FALSE,col.names=FALSE);

setwd("/home/limeng/splicingSNP/annotation/");

intronBedSort<-generatingIntronSplicingRegion("refgene.bed",2+1,1,2,1+1);
options("scipen"=1000, "digits"=4);
write.table(intronBedSort,file="splicingRegionBed3.bed",row.names=FALSE,quote=FALSE,col.names=FALSE);

