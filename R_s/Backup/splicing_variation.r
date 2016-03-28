library(GenomicRanges)
library(stringr)

vcfdata<-read.table("/home/limeng/introSNPdisease/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.clean.snp",sep="\t",header=FALSE,as.is=TRUE,nrow=100000)
colnames(vcfdata)<-c("chr","pos","ID","ref","alt","score","filter","annotation");


vcf <- with(vcfdata, GRanges(str_c("chr",chr), IRanges(pos, pos),NA, score,filter,annotation));
gtfdata<-read.table("/home/limeng/drugresistant/RNA-seq/mRNA/annotation/hg19refflat.gtf",sep="\t",as.is=TRUE,header=FALSE);
colnames(gtfdata)<-c("chr","source","region","beg","end","score","strand","frame","annotation");
gtf <- with(gtfdata, GRanges(chr, IRanges(beg, end),strand,score,annotation));

targetVCF <-  subsetByOverlaps(vcf,gtf);
vcfhits<-findOverlaps(vcf,gtf);
vcfgtf<-cbind(as.data.frame(vcf[queryHits(vcfhits)]),as.data.frame(gtf[subjectHits(vcfhits)]))




data<-read.table(file="/home/limeng/introSNPdisease/neutral.snp.ref.mut.seq.proximity",sep="\t",header=FALSE,as.is=FALSE,skip=1)
colnames(data)<-c("coordinate"  ,"transcript",	"snp",	"5'proximity",	"3'proximity", "maf",	"asn_maf",	"amr_maf",	"afr_maf",	"eur_maf")
dataExonSplicing<-data[(data[,4]==0)|(data[,5]==0),]

dataExonSplicingAF<-dataExonSplicing[(dataExonSplicing[,7]>=0.05)|(dataExonSplicing[,8]>=0.05)|(dataExonSplicing[,9]>=0.05)|(dataExonSplicing[,10]>=0.05),]
write.table(dataExonSplicingAF[,c(3,7,8,9,10)],file="/home/limeng/introSNPdisease/exonSNP",sep="\t",row.names=FALSE)


setwd("/home/limeng/splicingSNP/1000genomedata/")
dataIntron<-read.table("/home/limeng/splicingSNP/1000genomedata/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.clean.annovar.snp.variant_function_splicing~",header=FALSE,as.is=TRUE,sep="\t");
dataIntron[,8]<-str_c(dataIntron[,8],";")


result<-cbind(str_extract(dataIntron[,8],"ASN_AF=[^;]*.;"),str_extract(dataIntron[,8],"AMR_AF=[^;]*.;"),str_extract(dataIntron[,8],"AFR_AF=[^;]*.;"),str_extract(dataIntron[,8],"EUR_AF=[^;]*.;") )
#result<-apply(result,c(1,2),function(x){return( str_sub(x,8,(nchar(x)-1)) );     })


result<-apply(result,c(1,2),function(x){x=str_sub(x,8,(nchar(x)-1)); if(is.na(x)){ return(0)} else{return(as.numeric(x))} ;    }  )
resultAF<-(result[,1]>=0.005)|(result[,2]>=0.005)|(result[,3]>=0.005)|(result[,4]>=0.005) ;

dataIntronAFBig005<-dataIntron[resultAF,]
write.csv(dataIntronAFBig005,file="big005.csv",quote=FALSE)

resultAnno<-cbind(dataIntronAFBig005[,3],dataIntronAFBig005[,4],dataIntronAFBig005[,2],dataIntronAFBig005[,6],dataIntronAFBig005[,7] )
colnames(resultAnno)<-c("chr","position","id","ref","alt");

write.table(resultAnno,file="1000genomesplcingIntronAF",quot=FALSE,sep="\t",row.names=FALSE)


for(i in 0:((nrow(result)/40))  ){
    
  for(j in 1:40){
    write.table(resultAnno[(i*40+j):(i*40+j+39),],file=str_c("forPSI",i),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
    
    
  }
  
}



