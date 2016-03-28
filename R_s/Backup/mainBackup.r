
findHGMDExons<-function(exonPos){
  bed<-read.table("/home/limeng/splicingSNP/annotation/refgene.bed",sep="\t"   ,header=FALSE,as.is=TRUE);
  colnames(bed)<-c("chr","transcriptBeg","transcriptEnd","id","score","strand","thickStart","thickEnd","itemRgb","exonCount","blockSizes","blockStart")
  anno<-strsplit(exonPos,"\\.");
  
  transcriptID<-sapply(anno,"[",1);
  exonIndex<-as.numeric(sapply(anno,"[",2));
  
  
  exonsStart<-c();
  exonsEnd<-c();
  chrosome<-c();
  snpPos<-c();
  
  for(i in 1:length(transcriptID)){
    
    line<-which(bed[,"id"]==transcriptID[i]);
    if(length(line)==0){
      cat(str_c("can't find the current transcript ID:",transcriptID[i],"\n") );
      
      next;
    }
    
    transcript<-bed[line,];
    
    mTranscriptStart<-transcript[1,"transcriptBeg"];
    mTranscriptEnd<-transcript[1,"transcriptEnd"];
    mStrand<-transcript[1,"strand"];
    mExonCount<-as.numeric(transcript[1,"exonCount"]);
    mBlockSize<-as.numeric(strsplit(as.character(transcript[1,"blockSizes"]),",")[[1]] );
    mBlockStart<-as.numeric(strsplit(as.character(transcript[1,"blockStart"]),",")[[1]] );
    end<-mBlockStart+mBlockSize-1;
    
    mExonIndex<-exonIndex[i];
    if(mStrand=="-")
      mExonIndex<-mExonCount-mExonIndex+1;
    
    targetExonBeg<-mTranscriptStart+mBlockStart[mExonIndex]+1;
    targetExonEnd<-mTranscriptStart+end[mExonIndex]+1;
    
    
    if(length(targetExonBeg)==0||length(targetExonEnd)==0||is.null(targetExonEnd)||is.null(targetExonBeg)||is.na(targetExonBeg)||is.na(targetExonEnd) ){
      cat(str_c("can't find the exon:",transcriptID[i],"\n"))
      next;
      
    }
    
    exonsStart<-c(exonsStart,targetExonBeg);
    exonsEnd<-c(exonsEnd,targetExonEnd);
    chrosome<-c(chrosome,transcript[1,"chr"]);
    snpPos<-c(snpPos,)
    
  }
  
  
  xinjun<-str_c(chrosome,":",exonsStart,"-",exonsEnd);
  #not unique!! shouldn't be uniqued
  #xinjun<-unique(xinjun)
  
  write.table(xinjun,file="hgmdXinjunInput.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",");
  
  
}

