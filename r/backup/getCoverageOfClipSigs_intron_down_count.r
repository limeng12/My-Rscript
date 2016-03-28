setwd("E:\\limeng\\RBP_network\\r");
library(GenomicRanges);
library(dplyr)
library(reshape2)
library(stringr)

setwd("E:\\limeng\\RBP_network\\r")

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

exonId<-sapply(strsplit(dataAnno2[,"attr_miso"],";"),"[",1);

start_pos<-as.numeric(sapply(strsplit(exonId,":"),"[",6) );
end_pos<-as.numeric(sapply(strsplit(exonId,":"),"[",8) );

start_neg<-as.numeric(sapply(strsplit(exonId,":"),"[",6) );
end_neg<-as.numeric(sapply(strsplit(exonId,":"),"[",2) );

for(i in 1:length(dataAnno2[,"strand"]) ){
  if(dataAnno2[i,"strand"]=="+"){
    dataAnno2[i,"start"]<-start_pos[i];
    dataAnno2[i,"end"]<-end_pos[i];

  }else{
    dataAnno2[i,"start"]<-start_neg[i];
    dataAnno2[i,"end"]<-end_neg[i];
  }
}


dataAnno2[,"gene"]<-sapply(strsplit(dataAnno2[,"attr"],";"),
                          function(x){  return(paste0(x[grepl("gene_name",x)], x[grepl("exon_id",x)]) ) ;  }  );
dataAnno2<-dataAnno2[!duplicated(dataAnno2[,c("chr","start","end")]),];

                          
dataAnnoRange<-with(dataAnno2,
                    GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),
                            strand=strand,gene=gene,functions=functions,exonId=exonId));


files<-list.files("../clipbigBed/");
filesName<-sapply(strsplit(files,"\\."),"[",1);

clipRangesList<-GRangesList();
for(i in filesName){
  oneClip<-read.table(paste0("../clipbigBed/",i,".bigBed.bed"),header=F,as.is=T,sep="\t");
  colnames(oneClip)<-c("chr","start","end","gene","exon_number","strand","V7","V8","score","count");
  oneClip[,"start"]<-oneClip[,"start"]+1
  oneClipRange<-with(oneClip,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand));
  clipRangesList[[i]]<-oneClipRange;
}

#mtch<-findOverlaps(dataAnnoRange,clipRangesList);
geneClipMatch<-data.frame(genename=c(),clipname=c(),start=c(),end=c(),functions=c(),exonId=c(),stringsAsFactors = F);

for(i in 1:length(clipRangesList) ){
  oneMtch<-findOverlaps(dataAnnoRange,clipRangesList[[i]]);
  geneNames<-mcols(dataAnnoRange)[queryHits(oneMtch),"gene"];
  start<-start(dataAnnoRange)[queryHits(oneMtch)];
  end<-end(dataAnnoRange)[queryHits(oneMtch)];
  functions<-mcols(dataAnnoRange)[queryHits(oneMtch),"functions"];
  exonId<-mcols(dataAnnoRange)[queryHits(oneMtch),"exonId"];
  
  clipName<-names(clipRangesList)[i];
  geneClipMatch<-rbind(geneClipMatch,
                       data.frame(geneNames,rep(clipName,length(geneNames) ), start, end,functions,exonId) ); 
  #mcols(clipRangesList[[i]])[subjectHits(clipRangesList[[i]]),"gene"];
}

geneClipMatch[,1]<-as.character(geneClipMatch[,1]);
geneClipMatch[,2]<-as.character(geneClipMatch[,2]);

colnames(geneClipMatch)<-c("gene_name","clip_name","start","end","functions","exonId");


fileNameTargetMap<-read.table("../anno/clipNameTargetMap.txt",header=F,as.is=T,sep="\t");
colnames(fileNameTargetMap)<-c("target","fileName"); 
geneClipMatch2<-inner_join(geneClipMatch,fileNameTargetMap,by=c("clip_name"="fileName") );


b<-dcast(geneClipMatch2,gene_name+start+end+functions+exonId~target);
write.table(b,file="gene_clip_target.tsv",sep="\t");

d<-dcast(geneClipMatch2,gene_name+start+end+functions+exonId~target,
         function(x){if(length(x)!=0){ return(length(x))};return(0)});
write.table(d,file="gene_clip_target.tsv",sep="\t");

e_intron_down<-d[!is.na(d[,"functions"]),];
#e<-[,"exonId"]<-e[,"exonId"]
write.table(e_intron_down,file="gene_clip_target_se_intron_down.tsv",sep="\t");

