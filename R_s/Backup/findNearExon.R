
findTheExonsAndValidateTheGTAGConstraints<-function(snp){
  chrPosList<-strsplit(snp,":")
  chr<-(sapply(chrPosList,"[",1))
  pos<-as.numeric(sapply(chrPosList,"[",2) )
  vcf<-GRanges(chr, IRanges(pos, pos),"*")
  
  bedData<-read.table("/home/limeng/splicingSNP/annotation/refgene.bed",sep="\t",header=FALSE,as.is=TRUE);
  colnames(bedData)<-c("chr","start","end","id","score","strand","thickStart","thickEnd","color","blockCount","blockSizes","blockStart");
  bed <- with(bedData, GRanges(chr, IRanges(start, end),strand,blockCount,blockSizes,blockStart));
  
  
  vcfhits<-findOverlaps(vcf,bed);
  wrongAnnotationSNPs<-snp[-unique(queryHits(vcfhits))]
  cat(str_c("!!wrong annotation snps:",wrongAnnotationSNPs,"\n"))
  
  
  vcfFrame<-as(vcf[queryHits(vcfhits)],"data.frame");
  bedFrame<-as(bed[subjectHits(vcfhits)],"data.frame");
  colnames(vcfFrame)<-str_c("snp",colnames(vcfFrame));
  colnames(bedFrame)<-str_c("bed",colnames(bedFrame));
  
  vcfbed<-cbind(vcfFrame,bedFrame);
  
  
  bedblockEnds<-apply(vcfbed,1,function(x){
    size<-as.numeric(strsplit(x["bedblockSizes"],",")[[1]]); 
    start<-as.numeric(strsplit(x["bedblockStart"],",")[[1]]);
    end<-start+size-1;                           
    return(paste(end,",",collapse=""))
    
  });
  
  vcfbed<-cbind(vcfbed,bedblockEnds);
  
  refAlells<-as.character(getSeq(Hsapiens,vcfbed[,"snpseqnames"],vcfbed[,"snpstart"],vcfbed[,"snpstart"] ));
  vcfbed<-cbind(vcfbed,refAlells);
  
  
  result<-apply(vcfbed,1,function(x,desc){
    snppos<-as.numeric(x["snpstart"]);
    transcriptPos<-as.numeric(x["bedstart"]);
    transcriptEndPos<-as.numeric(x["bedend"]);
    
    
    #cat(str_c(class(snppos)," ",class(transcriptPos),"\n"))
    relativeSnpPos<-snppos-transcriptPos-1;#conver to 0-based coordinate
    #if(x["bedstrand"]=="-"){
    #  transcriptLength<-as.numeric(x["bedend"])-as.numeric(x["bedstart"])+1;
    #  relativeSnpPos<-transcriptLength-relativeSnpPos-1;  
    #}  
    
    
    exonStartPos<-as.numeric(strsplit(x["bedblockStart"],",")[[1]] );  
    exonEndPos<-as.numeric(strsplit(x["bedblockEnds"],",")[[1]] );
    #cat(exonStartPos,"#",exonEndPos,"\n")
    
    count<-as.numeric(x["bedblockCount"])
    #cat(str_c("count:",count))
    strand<-x["bedstrand"];
    
    #refAlell<-as.character(getSeq(Hsapiens,x["snpseqnames"],snppos,snppos));
    refAlell<-x["refAlells"];
    
    #cat(str_c(refAlell,"\n") );
    
    for(i in 1:count){
      disToStart<-relativeSnpPos-exonStartPos[i];
      disToEnd<-relativeSnpPos-exonEndPos[i];
      #cat(str_c(disToStart,"\n"))
      #cat(str_c(disToEnd,"\n"))
      
      
      if(abs(disToStart)==2||abs(disToStart)==1){
        distance<-abs(disToStart);
        #cat(str_c(" ",snpNumber))
        
        if(strand=="+"){
          if((distance==1)&&(refAlell!='G'))
            break;
          if((distance==2)&&(refAlell!='A'))
            break;
          
          return(str_c(x["snpseqnames"],",",x["snpstart"],",",refAlell,",","as",",",-abs(disToStart),",",strand,",",transcriptPos+exonStartPos[i]+1,",",transcriptPos+exonEndPos[i]+1));
          
        }
        
        if(abs(disToStart)=="-"){
          
          
          if((distance==1)&&(refAlell!='C'))
            break;
          if((distance==2)&&(refAlell!='A'))
            break;
          
          
          return(str_c(x["snpseqnames"],",",x["snpstart"],",",refAlell,",","ds",",",abs(disToStart),",",strand,",",transcriptPos+exonStartPos[i]+1,",",transcriptPos+exonEndPos[i]+1));
        }
      }
      
      
      if(abs(disToEnd)==2||abs(disToEnd)==1){
        #cat(str_c(" ",snpNumber))
        
        distance<-abs(disToEnd);
        
        if(strand=="+"){
          if((distance==1)&&(refAlell!='G'))
            break;
          if((distance==2)&&(refAlell!='T'))
            break;
          
          
          return(str_c(x["snpseqnames"],",",x["snpstart"],",",refAlell,",","ds",",",abs(disToEnd),",",strand,",",transcriptPos+exonStartPos[i]+1,",",transcriptPos+exonEndPos[i]+1));
        }
        
        if(strand=="-"){
          if((distance==1)&&(refAlell!='C'))
            break;
          if((distance==2)&&(refAlell!='T'))
            break;
          
          
          return(str_c(x["snpseqnames"],",",x["snpstart"],",",refAlell,",","as",",",-abs(disToEnd),",",strand,",",transcriptPos+exonStartPos[i]+1,",",transcriptPos+exonEndPos[i]+1));
        }
        
      }
      
      
    }
    
    return("");
    
    
  }  )
  
  
  result<-result[result!=""]
  #cat(result,file="tmpResult.txt",sep="\n");
  
  splitData<-strsplit(result,",");
  resultMatrix<-matrix(ncol=8,nrow=length(splitData));
  colnames(resultMatrix)<-c("chr","pos","ref","type","location","strand","exonBegPos","exonEndPos")
  
  for(i in 1:length(splitData)){
    resultMatrix[i,]<-splitData[[i]]
    
  }
  
  #sortByExonBegPos;
  #to make sure that HGMD and 1000 genome will pick the same exons;
  #resultMatrix<-resultMatrix[order(resultMatrix[,c("exonBegPos","exonEndPos") ]),]
  resultMatrix<-resultMatrix[order(resultMatrix[,"exonEndPos"],resultMatrix[,"exonBegPos"]),]
  
  #chrPos<-str_c(resultMatrix[,1],":",as.numeric(resultMatrix[,2]))
  #resultMatrixUnique<-resultMatrix[!duplicated(chrPos),]
  resultMatrixUnique<-resultMatrix;
  
  resultFormat<-str_c(resultMatrixUnique[,"chr"],":",as.numeric(resultMatrixUnique[,"pos"]),"\t",resultMatrixUnique[,"ref"],
                      "\t","N","\t",resultMatrixUnique[,"strand"],"\t",resultMatrixUnique[,"exonBegPos"],
                      "\t",resultMatrixUnique[,"exonEndPos"]);
  
  
  #cat(resultFormat,file="intronSnp.csv",sep="\n")
  
  
  
  #exonPos<-str_c(resultMatrixUnique[,"chr"],":",as.numeric(resultMatrixUnique[,"exonBegPos"]),"-",resultMatrixUnique[,"exonEndPos"])
  #resultMatrixUnique<-resultMatrixUnique[!duplicated(exonPos),];
  #write.csv(resultMatrixUnique,file="resultMatrixUnique.csv",row.names=FALSE)
  write.csv(resultMatrixUnique,file="resultMatrixUnique.csv",row.names=FALSE)
  
  
  exonCenterPos<-round(as.numeric(resultMatrixUnique[,"exonBegPos"])/2+as.numeric(resultMatrixUnique[,"exonEndPos"])/2)
  
  refs<-as.character(getSeq(Hsapiens,resultMatrixUnique[,"chr"],exonCenterPos,exonCenterPos));
  
  
  xinjunInput<-str_c(resultMatrixUnique[,"chr"],":",
                     exonCenterPos,"\t",
                     refs,"\t","T","\t",resultMatrixUnique[,"strand"],"\t",resultMatrixUnique[,"exonBegPos"],
                     "\t",resultMatrixUnique[,"exonEndPos"]
                     
  )
  
  cat(xinjunInput,file="xinjunIntronSnp.csv",sep="\n")
  
  xinEvolution<-str_c(resultMatrixUnique[,"chr"],":",resultMatrixUnique[,"exonBegPos"],"-",resultMatrixUnique[,"exonEndPos"])
  xinEvolution<-unique(xinEvolution)
  cat(xinEvolution,file="evolutionScore.csv",sep="\n")
  
}



findNearestExon<-function(snpFileName){
  snp<-read.table(snpFileName,header=TRUE,as.is=TRUE);
  chr<-snp[,1];
  
  #chr<-str_c("chr",snp[,1]);
  pos<-snp[,2];
  #chr<-(sapply(chrPosList,"[",1))
  #pos<-as.numeric(sapply(chrPosList,"[",2) )
  vcf<-GRanges(chr, IRanges(pos, pos),"*")
  
  bedData<-read.table("/home/limeng/splicingSNP/annotation/refgene.bed",sep="\t",header=FALSE,as.is=TRUE);
  colnames(bedData)<-c("chr","start","end","id","score","strand","thickStart","thickEnd","color","blockCount","blockSizes","blockStart");
  bed <- with(bedData, GRanges(chr, IRanges(start, end),strand,blockCount,blockSizes,blockStart));
  
  vcfhits<-findOverlaps(vcf,bed);
  wrongAnnotationSNPs<-snp[-unique(queryHits(vcfhits))]
  cat(str_c("!!wrong annotation snps:",wrongAnnotationSNPs,"\n"))
  
  
  vcfFrame<-as(vcf[queryHits(vcfhits)],"data.frame");
  bedFrame<-as(bed[subjectHits(vcfhits)],"data.frame");
  colnames(vcfFrame)<-str_c("snp",colnames(vcfFrame));
  colnames(bedFrame)<-str_c("bed",colnames(bedFrame));
  
  vcfbed<-cbind(vcfFrame,bedFrame);
  
  
  bedblockEnds<-apply(vcfbed,1,function(x){
    size<-as.numeric(strsplit(x["bedblockSizes"],",")[[1]]); 
    start<-as.numeric(strsplit(x["bedblockStart"],",")[[1]]);
    end<-start+size-1;                           
    return(paste(end,",",collapse=""))
    
  });
  
  vcfbed<-cbind(vcfbed,bedblockEnds);
  #cat(colnames(vcfbed)) 
  targetExons<-c();
  targetExonsAndSNPs<-matrix(ncol=4,nrow=0);
  
  #result<-apply(vcfbed,1,function(x,desc){
  for(i in 1:nrow(vcfbed)){
    x=vcfbed[i,]
    #cat(str_c("class of vcfbed line:",class(x)))
    
    chrosome=as.character(x[1,"bedseqnames"]);
    snppos<-as.numeric(x[1,"snpstart"]);
    transcriptPos<-as.numeric(x[1,"bedstart"]);
    transcriptEndPos<-as.numeric(x[1,"bedend"]);
    
    relativeSnpPos<-snppos-transcriptPos-1;#conver to 0-based coordinate
    
    #cat(str_c("class of bed:",class(x[1,"bedblockStart"]) ) )
    
    exonStartPoss<-as.numeric(strsplit(x[1,"bedblockStart"],",")[[1]] );  
    exonEndPoss<-as.numeric(strsplit(as.character(x[1,"bedblockEnds"]),",")[[1]] );
    #cat(exonStartPos,"#",exonEndPos,"\n")
    
    count<-as.numeric(x[1,"bedblockCount"])
    #cat(str_c("count:",count))
    strand<-x[1,"bedstrand"];
    
    exonIndex=-1;
    minToExonStart<-min(abs(relativeSnpPos-exonStartPoss) )
    minToExonEnd<-min(abs(relativeSnpPos-exonEndPoss) )
    
    if(minToExonStart>3&&minToExonEnd>3){
#cat(str_c("current transcript is not near the splicing site! distanceToExonStart:",minToExonStart,"distanceToExonsEnd:",minToExonEnd,"\n"  )  )
      next;
    }
      
    if(minToExonStart<minToExonEnd){
      exonIndex=which.min(abs(relativeSnpPos-exonStartPoss))
    }else{
      exonIndex=which.min(abs(relativeSnpPos-exonEndPoss))
    }
    
    exonBegPos<-exonStartPoss[exonIndex]+transcriptPos;
    exonEndPos<-exonEndPoss[exonIndex]+transcriptPos;
    
    #targetExons<-c(targetExons,str_c(chrosome,":",exonBegPos,"-",exonEndPos));
    targetExonsAndSNPs<-rbind(targetExonsAndSNPs,c(chrosome,snppos,exonBegPos,exonEndPos) ) 
    #cat(chrosome)
    #cat("\n")

  }
  
  write.csv(targetExonsAndSNPs,file="nearest exon with snps.csv",quote=FALSE,row.names=FALSE);
  colnames(targetExonsAndSNPs)<-c("chrosome","snppos","exonBegPos","exonEndPos");
  

  duplicatedLabel<-duplicated(str_c(targetExonsAndSNPs[,"chrosome"],targetExonsAndSNPs[,"snppos"]) );
  targetExonsAndSNPs<-targetExonsAndSNPs[!duplicatedLabel,];
  
targetExons<-str_c(targetExonsAndSNPs[,"chrosome"],":",targetExonsAndSNPs[,"exonBegPos"],"-",targetExonsAndSNPs[,"exonEndPos"]);
  #targetExons<-unique(targetExons)
  return(targetExons);
  
}


setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20")
targetExons<-findNearestExon("1000genome.phase3.splicingjunction20.SNPs_only.AF10%.recode.exon.only.recode.vcf");
cat(targetExons,file="exonsNearSNPS.csv",sep="\n")


