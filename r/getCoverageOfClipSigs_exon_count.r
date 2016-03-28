library(GenomicRanges);
library(dplyr)
library(reshape2)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(plyr)

setwd("/home/limeng/Projects/RBP_network/r")
#setwd("E:\\limeng\\RBP_network\\r")

dataAnno<-read.table("../anno/Homo_sapiens.GRCh37.75.gtf",header=F,as.is=T,sep="\t");
colnames(dataAnno)<-c("chr","function","class","start","end","score","strand","id","attr");
dataAnno[,"chr"]<-paste0("chr",dataAnno[,"chr"]);

#only gene
dataAnno<-dataAnno[dataAnno[,"class"]=="exon",];

skipExon<-read.table("../anno/genes.gff",header=F,as.is=T,sep="\t");
colnames(skipExon)<-c("chr_miso","functions","class","start","end","score","strand_miso","id","attr_miso"); 
exonId<-sapply(strsplit(skipExon[,"attr_miso"],";"),"[",1);
start<-sapply(strsplit(exonId,":"),"[",5);
end<-sapply(strsplit(exonId,":"),"[",6);
skipExon[,"start"]<-as.numeric(start);
skipExon[,"end"]<-as.numeric(end);
skipExon[,"exonId"]<-exonId;

dataAnno2<-inner_join(dataAnno,skipExon,by=c("start"="start","end"="end","strand"="strand_miso") );

dataAnno2[,"start"]<-dataAnno2[,"start"];
dataAnno2[,"end"]<-dataAnno2[,"end"];


dataAnno2[,"gene"]<-sapply(strsplit(dataAnno2[,"attr"],";"),
                          function(x){  return(paste0(x[grepl("gene_name",x)], x[grepl("exon_id",x)]) ) ;  }  );

dataAnno2<-dataAnno2[!duplicated(dataAnno2[,c("chr","start","end")]),];

dataAnnoRange<-with(dataAnno2,
                    GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                            strand=strand,gene=gene,functions=functions,exonId=exonId));

files<-list.files("../clipbed/");
filesName<-sapply(strsplit(files,"\\."),"[",1);

clipRangesList<-GRangesList();
for(i in filesName){
  oneClip<-read.table(paste0("../clipbed/",i,".bed"),header=TRUE,as.is=TRUE,sep="\t");
  colnames(oneClip)<-c("chr","start","end","gene","exon_number","strand","V7","V8","score","count");
  oneClip[,"start"]<-oneClip[,"start"]+1

  oneClipRange<-with(oneClip,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand));
  clipRangesList[[i]]<-oneClipRange;
}

#mtch<-findOverlaps(dataAnnoRange,clipRangesList);
geneClipMatch<-data.frame(genename=c(),clipname=c(),start=c(),end=c(),
   functions=c(),exonId=c(),start_rbp=c(),
   end_rbp=c(),chr_rbp=c()
   ,strand_rbp=c(),stringsAsFactors = F);
   
for(i in 1:length(clipRangesList) ){
  oneMtch<-findOverlaps(dataAnnoRange,clipRangesList[[i]]);
  geneNames<-mcols(dataAnnoRange)[queryHits(oneMtch),"gene"];
  start<-start(dataAnnoRange)[queryHits(oneMtch)];
  end<-end(dataAnnoRange)[queryHits(oneMtch)];
  functions<-mcols(dataAnnoRange)[queryHits(oneMtch),"functions"];
  exonId<-mcols(dataAnnoRange)[queryHits(oneMtch),"exonId"];
  
  start_rbp<-start(clipRangesList[[i]])[subjectHits(oneMtch)];
  end_rbp<-end(clipRangesList[[i]])[subjectHits(oneMtch)];
  chr_rbp<-seqnames(clipRangesList[[i]])[subjectHits(oneMtch)];
  strand_rbp<-strand(clipRangesList[[i]])[subjectHits(oneMtch)];
  
  clipName<-names(clipRangesList)[i];
  geneClipMatch<-rbind(geneClipMatch,
                       data.frame(geneNames,rep(clipName,length(geneNames) ), start, end,functions,exonId,start_rbp,end_rbp,chr_rbp,strand_rbp ) ); 
  #mcols(clipRangesList[[i]])[subjectHits(clipRangesList[[i]]),"gene"];
  
}

geneClipMatch[,1]<-as.character(geneClipMatch[,1]);
geneClipMatch[,2]<-as.character(geneClipMatch[,2]);

colnames(geneClipMatch)<-c("gene_name","target","start","end","functions","exonId","start_rbp","end_rbp","chr_rbp","strand_rbp");


#fileNameTargetMap<-read.table("../anno/clipNameTargetMap.txt",header=F,as.is=T,sep="\t");
#colnames(fileNameTargetMap)<-c("target","fileName"); 
#geneClipMatch2<-inner_join(geneClipMatch,fileNameTargetMap,by=c("clip_name"="fileName") );

if(FALSE){

ddply(geneClipMatch2,.(target),function(x){
  target_name<-x[1,"target"];
  target_name_file<-paste0("../rbp_bind_seq/",target_name,"_exon");
 
  cat("",file=target_name_file,append=F);
  label_contain<-c();
  
  for(i in 1:nrow(x)){
    if(is.element(x[i,"start_rbp"],label_contain) )
      next;
    label_contain<-c(label_contain,x[i,"start_rbp"])
    
    label<-paste0(">",x[i,"chr_rbp"],":",x[i,"start_rbp"],"-",x[i,"end_rbp"]);
    sequence<-getSeq(Hsapiens,x[i,"chr_rbp"],x[i,"start_rbp"],x[i,"end_rbp"]);
    if(nchar(sequence)<10)
      next;
    if(x[i,"strand_rbp"]=="-")
      sequence=reverseComplement(sequence);
    
    cat(paste0(label,"\n"),file=target_name_file,append=T);
    cat(paste0(sequence,"\n"),file=target_name_file,append=T)
    
  }
})

}

#b<-dcast(geneClipMatch2,gene_name+start+end+functions+exonId~target);
#write.table(b,file="gene_clip_target.tsv",sep="\t");

#coutn the number of target in each Event
d<-dcast(geneClipMatch,gene_name+start+end+functions+exonId~target,
         function(x){return(length(x)) });
write.table(d,file="gene_clip_target.tsv",sep="\t");

e_exon<-d[!is.na(d[,"functions"]),];
#e<-[,"exonId"]<-e[,"exonId"]
write.table(e_exon,file="gene_clip_target_se_exon.tsv",sep="\t");

