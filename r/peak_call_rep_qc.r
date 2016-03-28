library(GenomicRanges)
library(reshape2)
library(ggplot2)

#setwd("E:\\limeng\\RBP_network\\r");

RBP_sample_table<-read.table("../anno/clipNameTargetMap.txt",header=FALSE,sep="\t",as.is=TRUE);
colnames(RBP_sample_table)<-c("RBP", "sample");

RBP_names<-unique(RBP_sample_table[,"RBP"]);

RBP_sample_map<-list();
for(i in 1:length(RBP_names) ){
	One_RBP_samples<-RBP_sample_table[RBP_sample_table[,1]==RBP_names[i],2,drop=TRUE];
	RBP_sample_map[[RBP_names[i] ]]<-One_RBP_samples;
	
}

#RBPs_replicate_overlap<-data.frame(RBP_name<-c(),RBP_1_count<-c(),RBP_1_overLap_count<-c(),RBP_2_count<-c(),RBP_2_overLap_count<-c());
RBPs_replicate_overlap<-matrix(ncol=5,nrow=length(RBP_sample_map) );

for(i in 1:length(RBP_sample_map) ){
	RBP_name<-names(RBP_sample_map)[i];
	
	print(paste0(i," sample: ",RBP_name ) );
	
	sample_names<-RBP_sample_map[[i]];
	if(length(sample_names)<2)
		next;
	
	clip1<-read.table(paste0("../clipbigBed/",sample_names[1],".bigBed.bed"),header=FALSE,as.is=TRUE,sep="\t");
	numberOfClip1<-nrow(clip1);
	colnames(clip1)<-c("chr","start","end","gene","exon_number","strand","V7","V8","score","count");
	clip1[,"start"]<-clip1[,"start"]+1;
	clip1Range<-with(clip1,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand));	
	
	clip2<-read.table(paste0("../clipbigBed/",sample_names[2],".bigBed.bed"),header=FALSE,as.is=TRUE,sep="\t");
	numberOfClip2<-nrow(clip2);
	colnames(clip2)<-c("chr","start","end","gene","exon_number","strand","V7","V8","score","count");
	clip2[,"start"]<-clip2[,"start"]+1;
	clip2Range<-with(clip2,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end),strand=strand));	
	
	oneMtch<-findOverlaps(clip1Range,clip2Range);
	overLapNum1<-length(unique(queryHits(oneMtch)));
	overLapNum2<-length(unique(subjectHits(oneMtch)));
	
	clip1[,"start"]<-clip1[,"start"]-1;
	write.table(clip1[queryHits(oneMtch),], file=paste0("../clipbed/",RBP_name,".bed"),sep="\t",quote=FALSE,row.names=FALSE );
	
	#RBPs_replicate_overlap<-rbind(RBPs_replicate_overlap,c(RBP_name,numberOfClip1,numberOfClip2,overLapNum1,overLapNum2) );
	RBPs_replicate_overlap[i,]<-c(RBP_name,numberOfClip1,overLapNum1,numberOfClip2,overLapNum2)
}

colnames(RBPs_replicate_overlap)<-c("RBP","RBP_1_count","RBP_1_count_overlap","RBP_2_count","RBP_2_count_overlap");
RBPs_replicate_overlap_frame<-as.data.frame(RBPs_replicate_overlap);
#RBPs_replicate_overlap<-cbind(rownames(RBPs_replicate_overlap),RBPs_replicate_overlap);
RBPs_replicate_overlap_melt<-melt(RBPs_replicate_overlap_frame,id.vars="RBP");

RBPs_replicate_overlap_melt[,"value"]<-as.numeric(RBPs_replicate_overlap_melt[,"value"]);

pdf("clip_rep_qc.pdf",height=5,width=10);
p<-ggplot(RBPs_replicate_overlap_melt, aes(x=RBP, y=value,fill=variable)) + geom_bar(stat="identity",position="dodge")+
theme(axis.text.x = element_text(angle = 90, vjust = 1, 
    size = 6, hjust = 1),axis.text.y = element_text(angle = 0, vjust = 1, 
    size = 3, hjust = 1) );

print(p)
dev.off();


