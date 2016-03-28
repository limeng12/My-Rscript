library(GenomicRanges)
library(ShortRead)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)

getSeqFromUCSC<-function(snpPos){
  snpSplit<-strsplit(snpPos,":")
  
  
}


setwd("/home/limeng/introSNPdisease/introSNPdisease")
hgmdSNP<-read.table("Yunlong_HGMD_INTRONIC_SPLICING_MORT_14_2.csv",sep="\t",header=TRUE,as.is=TRUE,quote="",comment.char="")

hgmdSNPAs<-hgmdSNP[(hgmdSNP$type=="as")&(hgmdSNP$location==1),]
hgmdSNPDs<-hgmdSNP[(hgmdSNP$type=="ds")&(hgmdSNP$location==-1),]

hgmdExonSnp<-rbind(hgmdSNPAs,hgmdSNPDs);
hgmdExonSnpPos<-str_c(hgmdExonSnp$hg19_chromosome,":",hgmdExonSnp$hg19_coordinate)
cat(str_c("duplicate SNPs in HGMD:",hgmdExonSnpPos[duplicated(hgmdExonSnpPos)],"\n"))

hgmdExonSnpPos<-unique(hgmdExonSnpPos);


#intronSNP<-read.table("intronSNP",header=TRUE,sep="\t",as.is=TRUE)
one1000ExonSNP<-read.csv("exonSNP.csv",header=TRUE,sep=",",as.is=TRUE);
#colnames(intronSNP)<-colnames(exonSNP)
one1000Snp<-one1000ExonSNP[,1];
cat(str_c("duplicate SNPs in 1000genome:",one1000Snp[duplicated(one1000Snp)],"\n"))

one1000Snp<-unique(one1000Snp);



intersectionSNP<-intersect(one1000Snp,hgmdExonSnpPos);

snp<-setdiff(union(hgmdExonSnpPos,one1000Snp),intersect(one1000Snp,hgmdExonSnpPos));
cat(str_c("number of intersection SNPs:",length(intersectionSNP)," after  remove the intersection SNPs:",length(snp),"\n"));
cat(str_c("number of hgmd SNPs:",length(hgmdExonSnpPos)-length(intersectionSNP),"\n"))
cat(str_c("number of 1000genome SNPs:",length(one1000Snp)-length(intersectionSNP),"\n") )


#snps<-rbind(intronSNP,exonSNP);
chrPosList<-strsplit(snp,":")

chr<-(sapply(chrPosList,"[",1))

pos<-as.numeric(sapply(chrPosList,"[",2) )

vcf<-GRanges(chr, IRanges(pos, pos),"*")

gtfdata<-read.table("refgene.gtf",sep="\t",as.is=TRUE,header=FALSE)
colnames(gtfdata)<-c("chr","source","region","beg","end","score","strand","frame","annotation");
gtf <- with(gtfdata, GRanges(chr, IRanges(beg, end),strand,score,annotation));

vcfhits<-findOverlaps(vcf,gtf);
wrongAnnotationSNPs<-snp[-unique(queryHits(vcfhits))]
cat(str_c("!!wrong annotation snps:",wrongAnnotationSNPs,"\n"))


vcfFrame<-as.data.frame(vcf[queryHits(vcfhits)])
gtfFrame<-as.data.frame(gtf[subjectHits(vcfhits)])
colnames(vcfFrame)<-str_c("snp",colnames(vcfFrame))
colnames(gtfFrame)<-str_c("gtf",colnames(gtfFrame))


vcfgtf<-cbind(vcfFrame,gtfFrame);



vcfgtfRight<-vcfgtf[apply(vcfgtf,1,function(x){
  if((x["snpstart"]==x["gtfstart"])|(x["snpstart"]==x["gtfend"]))
    return(TRUE);
  
  return(FALSE);
  
}),];

snppos<-str_c(vcfgtfRight$snpseqnames,vcfgtfRight$snpstart);
vcfgtfAddPos<-cbind(vcfgtfRight,":",snppos);

targetVCF<-vcfgtfAddPos[!duplicated(vcfgtfAddPos$snppos),];


annotation<-targetVCF$gtfannotation;
annotationList<-strsplit(annotation,";");
annotationVec<-sapply(annotationList,"[",1)

targetVCFFormat<-with(targetVCF,data.frame(snpseqnames,snpstart,gtfstart,gtfend,gtfstrand,str_sub(annotationVec,9,nchar(annotationVec))))
colnames(targetVCFFormat)<-c("chr","snpPos","region beg","region end","strand","transcript_id")
result<-apply(targetVCFFormat,1,function(x){
  if(x["snpPos"]==x["region beg"]){
    if(x["strand"]=="+")
      return(str_c(x["chr"],",",x["snpPos"],",","as",",",+1,",",x["strand"],",",x["region beg"],",",x["region end"]))
    if(x["strand"]=="-")
      return(str_c(x["chr"],",",x["snpPos"],",","ds",",",-1,",",x["strand"],",",x["region beg"],",",x["region end"]))
  }
  
  if(x["snpPos"]==x["region end"]){
    if(x["strand"]=="+")
      return(str_c(x["chr"],",",x["snpPos"],",","ds",",",-1,",",x["strand"],",",x["region beg"],",",x["region end"]) );
    if(x["strand"]=="-")
      return(str_c(x["chr"],",",x["snpPos"],",","as",",",+1,",",x["strand"],",",x["region beg"],",",x["region end"]) );
  }
  
  return("");
  
})

#targetVCFFormat<-cbind(targetVCFFormat,pos);
wrongExonSplicingSNPs<-targetVCFFormat[nchar(result)==0,];
cat("wrong exon splicing SNPs:")
write.table(wrongExonSplicingSNPs,row.names=FALSE,col.names=FALSE);
colnames<-paste("chr","pos","type","location","strand","exonBegPos","exonEndPos",sep=",",collapse="")


cat(colnames,file="exonSnp.csv",sep="\n");
cat(result,file="exonSnp.csv",sep="\n",append=TRUE);

data<-read.csv(result,header=FALSE,as.is=TRUE)
#write.csv(targetVCFFormat,file="AllExonSplicingSnp.csv",row.names=FALSE);





