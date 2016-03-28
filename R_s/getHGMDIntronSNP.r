library(GenomicRanges)
library(stringr)
library(Biostrings)
library(seqinr)
library('BSgenome.Hsapiens.UCSC.hg19');


# For each snp, we get the nearest exon
# threadshold is the snp's pos relative to the splicing junction.
# threshold+3 is the snp's nearest distence to the exon(splicing junction).
getSplicingJunctionAndInExon25bp<-function(threshold){
  cat("get nearest exons of a HGMD snps and unique the exons!\n");
  #splicing Junction +- 25bp
  setwd(hgmdRawDataPath)
  hgmdSNP<-read.table("Yunlong_HGMD_INTRONIC_SPLICING_MORT_14_2.tsv",
                      sep="\t",header=TRUE,as.is=TRUE,quote="",comment.char="");
  
  hgmdSNPAs<-hgmdSNP[(hgmdSNP$type=="as")&((hgmdSNP$location<=threshold)&(hgmdSNP$location>=-2)),];
  hgmdSNPDs<-hgmdSNP[(hgmdSNP$type=="ds")&((hgmdSNP$location>=-1*threshold)&(hgmdSNP$location<=2)),];
  
  hgmdIntronSnp<-rbind(hgmdSNPAs,hgmdSNPDs);

  
  hgmdIntronSnpPos<-str_c(hgmdIntronSnp$hg19_chromosome,":"
                          ,hgmdIntronSnp$hg19_coordinate,":"
                          ,hgmdIntronSnp$sub);
  
  #cat(str_c("duplicate SNPs in HGMD:",hgmdIntronSnpPos[duplicated(hgmdIntronSnpPos)],"\n") );
  #cat(str_c("number of duplicate SNPs in HGMD:",nrow(hgmdIntronSnpPos) ) );
  
  #7885
  cat(str_c("Number of SNPs in -2->",threshold,"=",length(unique(hgmdIntronSnpPos) ),"\n") );
  
  setwd(hgmdSnvSpliceJunc25);
  write.table(hgmdIntronSnpPos,file="HGMDSplicingSNPS25.csv",sep="\t",row.names=FALSE,quote=FALSE);
  
  #1302 number of genes
  cat(str_c("number of genes:",length(unique(hgmdIntronSnp$GENE) ),"\n")  );
  
  b=new("BED");
  mBed<-readFromBedFile(b,refGeneBedFilePath);
  snp<-data.frame(as.character(hgmdIntronSnp$hg19_chromosome),
                  as.numeric(hgmdIntronSnp$hg19_coordinate),
                  hgmdIntronSnp$hg19_strand,stringsAsFactors =FALSE);
  
  hgmdIds<<-getGeneIdFromSnp(mBed,snp);
  
  hgmdExonRegion<-getNearestExons(mBed,snp,threshold+3);
  #remove the same exon
  hgmdExonRegion<-unique(hgmdExonRegion);
  
  chr<-hgmdExonRegion[,"chrosome"];
  exonBegPos<-as.numeric(hgmdExonRegion[,"exonBegPos"] );
  exonEndPos<-as.numeric(hgmdExonRegion[,"exonEndPos"] );
  strand<-hgmdExonRegion[,"strand"];
  
  exonCenterPos<-round( (exonBegPos+exonEndPos)/2  );
  refSeq<-as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,exonCenterPos,exonCenterPos) );
  
  exonRegions<-str_c(chr,":",exonBegPos,"-",exonEndPos);
  xinJunInput<-str_c(chr,":",exonCenterPos,"\t",refSeq,"\t","T","\t",strand);
  
  cat(exonRegions,file="hgmdNearestExonRegions25.csv",sep="\n");# for debuging
  cat(xinJunInput,file="hgmd",sep="\n");
  
  #3725 in GT-AG region

}

#threshold is the max distance to a exon.
#getHGMDExonsInSpliceJunction();

threshold=-1;
getSplicingJunctionAndInExon25bp(threshold);




getHgmdProteinseq<-function(){
  hgmdIds<-unique(hgmdIds);
  cat(hgmdIds,sep="\n",file="hgmdId.txt");
  
  x<-read.fasta("/home/limeng/splicingSNP/annotation/refgeneseqhg19.fa",
                as.string=TRUE,seqtype="AA");
  
  
  nmMappingFile<-read.table(file="/home/limeng/splicingSNP/annotation/nm_np_mapping_hgmd",
                            header=FALSE,sep="\t",as.is=TRUE);
  hgmdNpIds<-nmMappingFile[nmMappingFile[,1] %in% hgmdIds,2]
  #idNoVersion<-sapply( strsplit(as.character(id(x) ),"\\." ),"[",1 );
  #result<-x[which(idNoVersion %in% npId) ];
  #npProtein<-data.frame( np=sapply( strsplit(as.character(id(x)),"\\."),"[",1 ),
  #                       seq=sread(x),stringsAsFactors = FALSE );
  
  npProtein<-data.frame(
    np=sapply( strsplit(names(sapply(x,"[",1)),"\\.") ,"[",1), 
    seq=unname(sapply(x,"[",1)),stringsAsFactors=FALSE );
  
  cat(file="hgmdProtein.fasta",append=FALSE);
  
  line=0;
  for(i in 1:nrow(npProtein) ){
    
    currentId<-npProtein[i,1];
    
    if( currentId %in% hgmdNpIds){
      cat(paste(">",currentId,"\n",npProtein[i,2],"\n",sep=""),file="hgmdProtein.fasta",append=TRUE);
      
    }
    
  }
  
}


#getHgmdProteinseq();

