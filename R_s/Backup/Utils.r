# group by a field and 
# keep the highest record in each group
#
#limeng Feb17 2015

library(plyr)
library(stringr)
library(ggplot2)

sortAndUniquePSI<-function(fileName){
  data<-read.table(fileName,as.is=TRUE,header=FALSE,sep="\t");
  #,"log_reg_score","reg_score_percentile"
  
  colnames(data)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT")
  
  sortGroupDataDup<-ddply(data,.(id),function(x){  return(  x[order(abs(x[,"dPSI"]),decreasing=TRUE), ]    )   }  )
  
  write.csv(sortGroupDataDup,file="sort group data test.csv",row.names=FALSE)
  
  duplicate<-duplicated(sortGroupDataDup[,"id"])
  sortGroupData<-sortGroupDataDup[!duplicate,]
  sortGroupData<-sortGroupData[order(abs(sortGroupData[,"dPSI"])),]
  
  write.csv(sortGroupData,file="High dPSI unique.csv",row.names=FALSE)

}

setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20/PSIResultAF10/AF10");
sortAndUniquePSI("allPsiResult.txt");


getExonPositionForPSIFile<-function(psiFile,bedFileName){
  psiFrame<-read.csv(file=psiFile,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT")
  anno<-strsplit(psiFrame[,2],":")
  
  transcriptID=sapply(anno,"[",4);# trancript id 
  exonsIndex=sapply(anno,function(x){exon<-x[[5]];return(as.numeric(str_sub(exon,5,nchar(exon) ) ) )  } );# exon Index
  cat(str_c("the first transcript is:",transcriptID[1],"\n"))
  cat(str_c("the first transcript's exon is number:",exonsIndex[1],"\n"))
  
  transcriptExons<-data.frame(transcriptID,exonsIndex,psiFrame[,"dPSI"],stringsAsFactors=FALSE);
  colnames(transcriptExons)<-c("transcriptID","exonIndex",dPSI);
  
  highPSITranscripts<-transcriptExons[transcriptExons[,"dPSI"]>=5,  ]
  lowPSITranscripts<-transcriptExons[transcriptExons[,"dPSI"]<=5,  ]
  
  
  mBed<-new BED();
  mBed<-readFromBedFile(b,"/home/limeng/splicingSNP/annotation/refgene.bed");
  
  lowPSIExonRegions<-getExonsRegion(mBed,lowPSIExonRegions[,"transcriptID"],lowPSIExonRegions[,"exonIndex"]);
  highPSITranscripts<-getExonsRegion(mBed,highPSITranscripts[,"transcriptID"],highPSITranscripts[,"exonIndex"]);
  
  cat(lowPSIExonRegions,file="low psi exons.csv",sep="\n")
  cat(highPSITranscripts,file="high psi exons.csv",sep="\n")

  
}


setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20/PSIResultAF10/AF10");
getExonPositionForPSI("High dPSI unique.csv","/home/limeng/splicingSNP/annotation/refgene.bed");


getExonPositionForPSI<-function(psiFile,bedFileName){
  bed<-read.table(file=bedFileName,sep="\t",as.is=TRUE);
  colnames(bed)<-c("chr","transcriptBeg","transcriptEnd","id","score","strand","thickStart","thickEnd","itemRgb","exonCount","blockSizes","blockStart")
  
  psiFrame<-read.csv(file=psiFile,header=TRUE,as.is=TRUE);
  anno<-strsplit(psiFrame[,2],":")
  
  id<-psiFrame[,1]
  chr=sapply(anno,"[",1)
  strand=sapply(anno,"[",3);
  transcriptID=sapply(anno,"[",4);
  exonsIndex=sapply(anno,function(x){exon<-x[[5]];return(as.numeric(str_sub(exon,5,nchar(exon) ) ) )  } );
  
  exonsStart<-c();
  exonsEnd<-c();
  chrosome<-c();
  psi<-c();
  psiTranscript<-c();
  psiID<-c();
  
  for(i in 1:nrow(psiFrame)){
    
    line<-which(bed[,"id"]==transcriptID[i]);
    if(length(line)==0){
      cat(str_c("can't find the current transcript ID:",transcriptID[i]," rs:",id[i],"\n") );
      
      next;
    }
    
    transcript<-bed[line,];
    
    mTranscriptStart<-transcript[1,"transcriptBeg"];
    mTranscriptEnd<-transcript[1,"transcriptEnd"];
    mStrand<-transcript[1,"strand"];
    mExonCount<-as.numeric(transcript[1,"exonCount"]);
    mBlockSize<-as.numeric(strsplit(as.character(transcript[1,"blockSizes"]),",")[[1]] ) ;
    mBlockStart<-as.numeric(strsplit(as.character(transcript[1,"blockStart"]),",")[[1]] ) ;
    end<-mBlockStart+mBlockSize-1;
    
    mExonIndex<-exonsIndex[i];
    if(mStrand=="-")
      mExonIndex<-mExonCount-mExonIndex+1;
    
    
    #cat(mBlockStart);cat(mExonIndex);
    targetExonBeg<-mTranscriptStart+mBlockStart[mExonIndex]+1;
    targetExonEnd<-mTranscriptStart+end[mExonIndex]+1;
    
    if(is.na(targetExonBeg)||is.na(targetExonEnd)){
      cat(str_c("can't find the exon:",transcriptID[i]," rs:",id[i])   )
      next;
    }
    
    exonsStart<-c(exonsStart,targetExonBeg);
    exonsEnd<-c(exonsEnd,targetExonEnd);
    chrosome<-c(chrosome,chr[i])
    psi<-c(psi,psiFrame[i,3])
    psiTranscript<-c(psiTranscript,psiFrame[i,2]);
    psiID<-c(psiID,id[i]);
    
  }
  
  exonsPsi<-data.frame(psi,chrosome,exonsStart,exonsEnd,stringAsFactors=FALSE);
  write.csv(cbind(exonsPsi,psiTranscript,psiID),file="psi with exon start and end.csv",row.names=FALSE,quote=FALSE);
  
  
  #for xinjun's program
  #exonsPsi<-cbind(psi,chrosome,exonsStart,exonsEnd);
  #write.csv(exonsPsi,file="psi with exon start and end.csv",row.names=FALSE,quote=FALSE)
  lowPsiExons<-exonsPsi[which(abs(exonsPsi[,1])<5),]
  highPsiExons<-exonsPsi[which(abs(exonsPsi[,1])>=5),]
  write.table(lowPsiExons,file="xinjunLowPsiExons.csv",quote=FALSE,row.names=FALSE)
  write.table(highPsiExons,file="xinjunHighPsiExons.csv",quote=FALSE,row.names=FALSE)
  
  
  xinjunLowExons<-with(lowPsiExons,str_c(chrosome,":",exonsStart,"-",exonsEnd)  );
  xinjunHighExons<-with(highPsiExons,str_c(chrosome,":",exonsStart,"-",exonsEnd)  );
  #write.csv(xinjun,file="xinJunProgramInput.csv",row.names=FALSE,quote=FALSE)
  #write.csv(xinjun,file="xinJunProgramInput.csv",row.names=FALSE,quote=FALSE)
  #png("psi.png")
  #hist
  #dev.off();
  
  xinjun<-str_c(chrosome,":",exonsStart,"-",exonsEnd);
  write.table(xinjun,file="xinJunProgramInput.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",");
  
  #lowPsiXinjun<-with(xinjunLowExons,str_c(chrosome,":",exonsStart,"-",exonsEnd) );
  write.table(xinjunLowExons,file="xinjunLowExons.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",");
  
  #highPsiXinjun<-with(xinjunHighExons,str_c(chrosome,":",exonsStart,"-",exonsEnd) );
  write.table(xinjunHighExons,file="xinjunHighExons.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",");
  
}

setwd("/home/limeng/splicingSNP/1000genome/SplicingJunction20/PSIResultAF10/AF10");
getExonPositionForPSI("High dPSI unique.csv","/home/limeng/splicingSNP/annotation/refgene.bed");


  
setwd("/home/limeng/splicingSNP/1000genomedata/PSI");
plotEvolution("lowPSI.txt.phylop","highPSI.txt.phylop")

getSNPFromPSI<-function(psiFile,inputToPSIFile){
  psiFrame<-read.csv(psiFile,header=TRUE,as.is=TRUE);
  
  id<-psiFrame[,1];
  anno<-strsplit(psiFrame[,2],":");
  chr=sapply(anno,"[",1);
  strand=sapply(anno,"[",3);
  transcriptID=sapply(anno,"[",4);
  
  inputPsi<-read.table(inputToPSIFile,sep="\t",header=FALSE,as.is=TRUE);
  inputPsi<-inputPsi[inputPsi[,3]!=".",]
  #inputPsiSelect<-inputPsi[(),];
  rownames(inputPsi)<-inputPsi[,3];
  inputPsiOrder<-inputPsi[id,];
  
  xinJunResult=str_c(chr,":",inputPsiOrder[,2],"\t",inputPsiOrder[,4],"\t",inputPsiOrder[,5],"\t",strand);
  write.csv(xinJunResult,file="xinJunNormalInput.csv",row.names=FALSE,quote=FALSE);
  
}

setwd("/home/limeng/splicingSNP/1000genomedata/PSIResult/AF10/");
getSNPFromPSI("High dPSI unique.csv","1000genome.phase3.splicingjunction20.SNPs_only.AF10%.recode.forPSI");
#psiFile<-"High dPSI unique.csv"
#inputToPSIFile<-"1000genome.phase3.splicingjunction20.SNPs_only.AF10%.recode.forPSI"


