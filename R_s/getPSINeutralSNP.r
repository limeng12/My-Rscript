# group by a field and 
# keep the highest record in each group
#
#limeng Feb17 2015
setwd(RPath);
library(plyr);
library(stringr);
library(ggplot2);
library('BSgenome.Hsapiens.UCSC.hg19');
source("BED.r");

#sort and unique a PSI result file, keep the highest |dPSI| for each snp.
sortAndUniquePSI<-function(fileName){
  cat("sort and unique a PSI result file, keep the highest |dPSI| for each snp\n");
  
  data<-read.table(fileName,as.is=TRUE,header=FALSE,sep="\t");
  #,"log_reg_score","reg_score_percentile"
  
  colnames(data)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
  
  sortGroupDataDup<-ddply(data,.(id),function(x)
    {  return(  x[order(abs(x[,"dPSI"]),decreasing=TRUE), ] ) } );
  
  write.csv(sortGroupDataDup,file="sort group data test.csv",row.names=FALSE);
  
  duplicate<-duplicated(sortGroupDataDup[,"id"])
  sortGroupData<-sortGroupDataDup[!duplicate,]
  sortGroupData<-sortGroupData[order(abs(sortGroupData[,"dPSI"])),];
  sortedFileName="High dPSI unique.csv"
  write.csv(sortGroupData,sortedFileName,row.names=FALSE);
  
  return(sortedFileName);
  
}


getExonPositionForPSIFile<-function(psiFile,bedFileName,threshold){
  psiFrame<-read.csv(file=psiFile,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
  anno<-strsplit(psiFrame[,"transcript"],":");
  
  transcriptID=sapply(anno,"[",4);# trancript id 
  exonsIndex=sapply(anno,function(x){exon<-x[[5]];
                                     return(as.numeric(str_sub(exon,5,nchar(exon) ) ) );  } );# exon Index
  
  cat(str_c("the first transcript is:",transcriptID[1],"\n"));
  cat(str_c("the first transcript's exon is number:",exonsIndex[1],"\n"));
  
  transcriptExons<-data.frame(transcriptID,exonsIndex,psiFrame[,"dPSI"],stringsAsFactors=FALSE);
  colnames(transcriptExons)<-c("transcriptID","exonIndex","dPSI");
  cat(str_c("number of |dPSI|>=10:",sum(abs(transcriptExons[,"dPSI"])>=10)  ) )
  
  highPSITranscripts<-transcriptExons[abs(transcriptExons[,"dPSI"])>=threshold,  ];
  lowPSITranscripts<-transcriptExons[abs(transcriptExons[,"dPSI"])<threshold,  ];
  a<-qplot(x=abs(highPSITranscripts[,"dPSI"]),geom="histogram",main="High dPSI histogram",xlab="|dpSI|");
  print(a)
  
  mBed<-new("BED");
  mBed<-readFromBedFile(b,refGeneBedFilePath);
  
  lowPSIExonRegions<-getExonsRegion(mBed,
                                    lowPSITranscripts[,"transcriptID"],
                                    lowPSITranscripts[,"exonIndex"]);
  
  highPSIExonRegions<-getExonsRegion(mBed,
                                     highPSITranscripts[,"transcriptID"],
                                     highPSITranscripts[,"exonIndex"]);
  
  
  cat(lowPSIExonRegions,file=str_c(threshold,"low.psi.exons.csv"),sep="\n");
  cat(highPSIExonRegions,file=str_c(threshold,"high.psi.exons.csv"),sep="\n");
  
  #ExonCentrenRegion
}


getNeutralExonCentralPos<-function(psiFileName,threshold){
  cat("get the center position of a exon using the PSI file, note this will unique the exons!\n");
  
  psiFrame<-read.csv(file=psiFileName,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
    
  #  |dPSI|>5
  psiFrame<-psiFrame[abs(psiFrame[,"dPSI"])>threshold,];
  
  anno<-strsplit(psiFrame[,"transcript"],":");
  transcriptID=sapply(anno,"[",4);# trancript id 
  exonsIndex=sapply(
    anno,function(x){exon<-x[[5]];
                     return(as.numeric(str_sub(exon,5,nchar(exon)  ) ) );  } );# exon Index
  
  strand=sapply(anno,"[",3);
  psi<-psiFrame[,"dPSI"];
  
  mBed=new("BED");
  
  mBed<-readFromBedFile(mBed,"/Users/mengli/Documents/projects/splicingSNP/annotation/refgene.bed");
  
  exonRegions<-getExonsRegion(mBed,transcriptID,exonsIndex);
  
  
  rsAF<-read.csv(file="rsAF.csv",header=FALSE,as.is=TRUE);
  rsAF<-rsAF[rsAF[,1]!=".",];
  row.names(rsAF)<-rsAF[,1];
  af<-rsAF[psiFrame[,"id"],2];
  #af<-af[!is.na(af) ]
  
  
  regions<-strsplit(exonRegions,":");
  chr_non_unique=sapply(regions,"[",1);
  beg_non_unique=as.numeric( sapply(regions,"[",2) );
  end_non_unique=as.numeric( sapply(regions,"[",3) );
  strand_non_unique=sapply(regions,"[",4);
  exonCenter_non_unique<-round( (beg_non_unique+end_non_unique)/2 );
  mappingFile<-data.frame(exonCenter=str_c(chr_non_unique,":",exonCenter_non_unique),
                          af=af,psiFrame[,"dPSI"],id=psiFrame[,"id"]  );
  
  mappingFile<-mappingFile[!duplicated(mappingFile[,1]),];
  colnames(mappingFile)<-c("id","AF","PSI");
  
  write.csv(mappingFile,file="mappingFile.csv",quote=FALSE,row.names=FALSE);
  
  #unique the exons
  exonRegions<-exonRegions[!(nchar(exonRegions)==0) ];
  regions<-strsplit(exonRegions,":");
  regions<-unique(regions);
  
  
  chr=sapply(regions,"[",1);
  beg=as.numeric( sapply(regions,"[",2) );
  end=as.numeric( sapply(regions,"[",3) );
  strand=sapply(regions,"[",4);
  exonCenter<-round( (beg+end)/2 );
  
  
  refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenter,exonCenter));
  
  result<-str_c(chr,":",exonCenter,"\t",refSeq,"\t","T","\t",strand,"\t",beg,"\t",end);
  
  #result<-str_c(exonCenterPos,targetSNPs[,"ref"],targetSNPs[,"alt"],,sep="\t");
  
  cat(result,file="normal",sep="\n");
  
}


getYaoqiNeutral<-function(){
  
  setwd(yaoNeutralPath);
  
  data<-read.table("mut_all_raw.txt",header=FALSE,as.is=TRUE,skip=2207);
  
  cat(str_c("number of neutral exons of yaoqi:",nrow(data),"\n"));
  
  snpSplit<-strsplit(data[,1],"\\|");
  chr<-sapply(snpSplit,"[",1);
  pos<-as.numeric(sapply(snpSplit,"[",2));
  
  snp<-data.frame(chr=chr,pos=pos,strand="*",stringsAsFactors=FALSE);
  #snp<-snp[snp[,"flag"]==0,]
  
  b=new("BED");
  mBed<-readFromBedFile(b,refGeneBedFilePath);
  yaoqiNearestExons<-getNearestExons(mBed,snp,3);
  write.csv(yaoqiNearestExons,"neutral exons near yaoqi.csv",row.names=FALSE,quote=FALSE);
  
  regions<-unique(yaoqiNearestExons);
  
  chr=regions[,1];
  beg=as.numeric( regions[,2] );
  end=as.numeric( regions[,3] );
  strand=regions[,4];
  exonCenter<-round( (beg+end)/2 );
  
  refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenter,exonCenter));
  
  result<-str_c(chr,":",exonCenter,"\t",refSeq,"\t","T","\t",strand,"\t",beg,"\t",end);
  
  #result<-str_c(exonCenterPos,targetSNPs[,"ref"],targetSNPs[,"alt"],,sep="\t");
  
  cat(result,file="normal",sep="\n");
  
}


getNeutralExonCentralPos2<-function(psiFileName,threshold,outputFileName){
  cat("get the center position of a exon using the PSI file, note this will unique the exons!\n");
  
  psiFrame<-read.csv(file=psiFileName,header=TRUE,as.is=TRUE);
  colnames(psiFrame)<-c("id","transcript","dPSI","dPSI_percentile","PSI_WT");
  #colnames(psiFrame)[10]<-"dPSI"; 
  #colnames(psiFrame)[14]<-"transcript"; 
  #colnames(psiFrame)[15]<-"exonIndex"; 
  
  #  |dPSI|>5
  psiFrame<-psiFrame[abs(psiFrame[,"dPSI"])>threshold,];
  
  anno<-strsplit(psiFrame[,"transcript"],":");
  #transcriptID<-psiFrame[,"transcript"];# trancript id 
  #exonsIndex<-psiFrame[,"exonIndex"];
  
  transcriptID<-sapply(anno,"[",4); 
  
  exonsIndex=sapply(
    anno,function(x){exon<-x[[5]];
                     return(as.numeric(str_sub(exon,5,nchar(exon)  ) ) );  } );# exon Index
  
  #strand=sapply(anno,"[",3);
  #psi<-psiFrame[,"dPSI"];
  
  mBed=new("BED");
  
  mBed<-readFromBedFile(mBed,"/Users/mengli/Documents/projects/splicingSNP/annotation/refgene.bed");
  exonRegions<-getExonsRegion(mBed,transcriptID,exonsIndex);
  
  #unique the exons
  exonRegions<-exonRegions[!(nchar(exonRegions)==0) ];
  regions<-strsplit(exonRegions,":");
  regions<-unique(regions);
  
  chr=sapply(regions,"[",1);
  beg=as.numeric( sapply(regions,"[",2) );
  end=as.numeric( sapply(regions,"[",3) );
  #strand=sapply(regions,"[",4);
  exonCenter<-round( (beg+end)/2 );
  
  #refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenter,exonCenter));
  #result<-str_c(chr,":",exonCenter,"\t",refSeq,"\t","T","\t",strand,"\t",beg,"\t",end);
  result<-str_c(chr,":",exonCenter); 
  
  cat(result,file=outputFileName,sep="\n");
}
